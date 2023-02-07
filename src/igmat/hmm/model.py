import os
import sys
import json
import math
import tempfile
import shutil
import traceback
from subprocess import Popen, PIPE

from igmat.alphabet import Alphabet
from igmat.hmm.result import Result
from igmat.hmm.mapper import Mapper

# Import the HMMER parser from the distributed version of Biopython.
try: 
  from Bio.SearchIO.HmmerIO import Hmmer3TextParser as HMMERParser
  from Bio.SearchIO.HmmerIO import Hmmer3TextIndexer
except ImportError as e:
  print('Unable to import Biopython. Please install Biopython module.')
  sys.exit(1)

def hmmer_path(name, hmmerpath=None):
  if hmmerpath:
    return os.path.join(hmmerpath, name)

  # Try to fetch path
  result = shutil.which(name)
  result = name if not result else result
  
  return result

class Model():

  def __init__(self, path, model='IMGT', hmmerpath=None):

    self.modelname = model
    self.modelpath = os.path.join(path, model)
    self.hmmerpath = hmmerpath
    self._dataset = {}

    if not self.__load__(path, model):
      raise Exception('Unable to correctly load HMMR model')

  def __load__(self, path, model):

    # Load the model information and store them in this class
    try:

      # Check if model file exists
      hmmerpath = os.path.join(path, model + '.hmm')
      if not os.path.exists(hmmerpath):
        raise Exception('Unable to find hmmer model.')

      # Run hmmer as a subprocess
      hmmstat = hmmer_path('hmmstat', self.hmmerpath)
      # hmmstat = "hmmstat"
      # if self.hmmerpath:
      #   hmmstat = os.path.join(self.hmmerpath, hmmstat)

      command = [ hmmstat, hmmerpath]
      process = Popen( command, stdout=PIPE, stderr=PIPE)

      # Parse the hmmstat output
      while True:
        line = process.stdout.readline().decode("utf-8")
        if not line:
          break

        if line[0] == '#':
          continue

        line = line.strip().split()
        if len(line) < 5:
          continue

        name = line[1]
        # species, chain = name.split('_')
        speciesList = name.split('_')
        chain = speciesList.pop()
        species = '_'.join(speciesList)
        self._dataset[name] = {
          'nseq': int(line[3]),
          'species': species,
          'chain': chain,
          'relent': float(line[6]), # Relative entropy
          'size': int(line[5])
        }

      # Load dataset info from the json file
      jsonpath = os.path.join(path, model + '.json')
      if not os.path.exists(jsonpath):
        raise Exception('Unable to find dataset json file')

      # Read alphabet data
      with open(jsonpath, 'r') as handle:
        data = json.load(handle)
        self.alpha = data['alphabet']

    except Exception as e:
        print('Error while opening HMMR model: %s' % e)
        return False

    return True

  def run(self, sequence, threshold=0, restrict=[]):

    # Run hmmer as a subprocess
    hmmscan = hmmer_path('hmmscan', self.hmmerpath)
    # hmmscan = "hmmscan"
    # if self.hmmerpath:
    #   hmmscan = os.path.join(self.hmmerpath, hmmscan)

    try:

      # Create a temporary directory
      dirpath = tempfile.mkdtemp()

      # Create a fasta file for all the sequences. Label them with their sequence index
      fasta_filename = os.path.join(dirpath, 'input.fasta')
      output_filename = os.path.join(dirpath, 'output.txt')
      with open(fasta_filename,'w') as handle:
        handle.write(">{name}\n{sequence}\n".format(
          name=sequence.getName(),
          sequence=Alphabet.reduce(sequence.getSequence(), self.alpha)
        ))

      command = [ hmmscan, "--domT", "0", "-o", output_filename, self.modelpath + '.hmm', fasta_filename]
      process = Popen(command, stdout=PIPE, stderr=PIPE  )
      _, pr_stderr = process.communicate()
      if pr_stderr:
          raise Exception(pr_stderr)

      # Parse the result
      parser = Hmmer3TextIndexer(output_filename)
      query = parser.get(0)

      domainMap = {}
      domainMask = ['*'] * sequence.getSize()

      # Iterate over the matches of the domains in order of their e-value (most significant first)
      for hsp in sorted(query.hsps, key=lambda x: x.evalue):

        # Skipping matches different to the top one
        nameList = hsp.hit_id.split('_')
        ctype = nameList.pop()
        species = '_'.join(nameList)

        # Check if the hit is an extension of a seed
        if hsp.hit_id in domainMap:
          skip_hit = False
          for hit_id in domainMap:
            if hit_id == hsp.hit_id:
              continue 

            if (hsp.query_start > domainMap[hit_id]['start'] and hsp.query_start < domainMap[hit_id]['end']) or (hsp.query_end > domainMap[hit_id]['start'] and hsp.query_end < domainMap[hit_id]['end']):
              skip_hit = True
              break

          if skip_hit:
            continue



        hitMask = [ '' if domainMask[i] == '*' else '*' for i in range(hsp.query_start, hsp.query_end)]
        overlaps = True if len(''.join(hitMask)) > 0 else False
        # if (hsp.hit_id not in domainMap) and overlaps:
        if overlaps:
          # print('Overlapping')
          continue

        if hsp.hit_id not in domainMap:
          domainMap[hsp.hit_id] = {
            'type': ctype,
            'start': -1,
            'end': -1,
            'evalue': None,
            'species': species,
            'hits': [],
            'list': []
          }

        # Check if in the exclude list
        if restrict and ctype not in restrict:
          logging.info('skipping chain {chain}'.format(chain=ctype))
          continue

        # This is a valid domain
        domainMask = [ ('*' if domainMask[i] == '*' and (i < hsp.query_start or i >= hsp.query_end) else '-') for i in range(sequence.getSize()) ]
        domainMap[hsp.hit_id]['start'] = min(hsp.query_start, domainMap[hsp.hit_id]['start']) if domainMap[hsp.hit_id]['start'] >= 0 else hsp.query_start
        domainMap[hsp.hit_id]['end'] = max(hsp.query_end, domainMap[hsp.hit_id]['end']) if domainMap[hsp.hit_id]['end'] >= 0 else hsp.query_end
        domainMap[hsp.hit_id]['evalue'] = min(hsp.evalue, domainMap[hsp.hit_id]['evalue']) if domainMap[hsp.hit_id]['evalue'] else hsp.evalue
        domainMap[hsp.hit_id]['hits'].append({
          'id': hsp.hit_id,
          'description': hsp.hit_description,
          'evalue': hsp.evalue,
          'bitscore': hsp.bitscore,
          'bias': hsp.bias,
          'start': hsp.query_start,
          'end': hsp.query_end
        })

        domainMap[hsp.hit_id]['list'].append(hsp)

      # No valid match found
      if not domainMap:
        return None

      resultList = []
      for key in domainMap:

        # Check the match seed evalue
        if domainMap[key]['evalue'] > 1e-10:
          continue
        
        # Generate a consensus of the found domains
        consensus = self.__hmm_align__(domainMap[key]['list'], query.seq_len, sequence)
        if not consensus:
          return None

        # Validate the alignment
        self.__hmm_validate__(consensus['state'], sequence)

        # Annotate the result
        alignment, annotation = self.__hmm_annotate__(consensus['state'], sequence)

        match = Result(alignment, consensus['start'], consensus['end'], annotation, domainMap[key]['type'], domainMap[key]['hits'])
        resultList.append(match)

      return resultList

    finally:
      shutil.rmtree(dirpath)
        
    # This shouldn't happen
    return None

  def __hmm_align__(self, hspList, seq_length, sequence):

    regionMap = {
      'FR1': {'start': 0, 'stop': 26},
      'CDR1': {'start': 26, 'stop': 38},
      'FR2': {'start': 38, 'stop': 55},
      'CDR2': {'start': 55, 'stop': 65},
      'FR3': {'start': 65, 'stop': 104},
      'CDR3': {'start': 104, 'stop': 117},
      'FR4': {'start': 117, 'stop': 127}
    }

    stateList = []
    for index, hsp in enumerate(hspList):

      # Extract the strings for the reference states and the posterior probability strings     
      reference_string = hsp.aln_annotation["RF"]
      state_string = hsp.aln_annotation["PP"]

      # print(hsp)
      # print(reference_string)
      # print(state_string)
      # print(hsp.hit_start, hsp.hit_end, hsp.query_start, hsp.query_end)

      # Extract the start an end points of the hmm states and the sequence
      # These are python indices i.e list[ start:end ] and therefore start will be one less than in the text file
      _hmm_start = hsp.hit_start
      _hmm_end = hsp.hit_end
       
      _seq_start = hsp.query_start
      _seq_end = hsp.query_end

      # Extact the full length of the HMM hit
      # species, ctype = hsp.hit_id.split('_')


      nameList = hsp.hit_id.split('_')
      ctype = nameList.pop()
      species = '_'.join(nameList)
      _hmm_length = self.size(species, ctype)

      # Handle cases where there are n terminal modifications.
      # In most cases the user is going to want these included in the numbered domain even though they are not 'antibody like' and 
      # not matched to the germline. Only allow up to a maximum of 5 unmatched states at the start of the domain
      # Adds a bug here if there is a very short linker between a scfv domains with a modified n-term second domain
      # Thus this is only done for the first identified domain ( hence order attribute on hsp )
      if index == 0 and _hmm_start and _hmm_start < 5: 
        n_extend = _hmm_start 
        if _hmm_start > _seq_start:
            n_extend = min( _seq_start , _hmm_start - _seq_start )
        state_string = '8'*n_extend + state_string  
        reference_string = 'x'*n_extend + reference_string
        _seq_start = _seq_start - n_extend
        _hmm_start = _hmm_start - n_extend

      # Handle cases where the alignment should be extended to the end of the j-element
      # This occurs when there a c-terminal modifications of the variable domain that are significantly different to germline
      # Extension is only made when half of framework 4 has been recognised and there is only one domain recognised.
      if len(hspList) ==1 and _seq_end < seq_length and (123 < _hmm_end < _hmm_length): # Extend forwards
        n_extend = min( _hmm_length - _hmm_end, seq_length - _seq_end )
        state_string = state_string + '8'*n_extend
        reference_string = reference_string + 'x'*n_extend
        _seq_end = _seq_end + n_extend
        _hmm_end = _hmm_end + n_extend
                      
      # Generate lists for the states and the sequence indices that are included in this alignment
      hmm_states = list(range(_hmm_start+1, _hmm_end+1))
      sequence_indices = list(range(_seq_start,  _seq_end))
      h, s = 0, 0 # initialise the current index in the hmm and the sequence
      
      state = []
      # iterate over the state string (or the reference string)
      for i in range( len(state_string) ):
        if reference_string[i] == "x": # match state
          state_type = "m"
        else: # insert state
          state_type = "i"
        
        if state_string[i] == ".": # overloading if deleted relative to reference. delete_state
          state_type = "d"
          sequence_index = None
          score = 0
        else:
          sequence_index = sequence_indices[s]
          score = (10 if state_string[i] == '*' else int(state_string[i]))

        # Store the alignment as the state identifier (uncorrected IMGT annotation) and the index of the sequence
        state.append({
          'hmm': hmm_states[h],
          'type': state_type,
          'idx': sequence_index,
          'score': score
        })

        # Updates to the indices         
        if state_type == "m":
          h+=1
          s+=1
        elif state_type == "i":
          s+=1
        else: # delete state
          h+=1

      # if len(stateList) == 0:
      stateList.append(state)

    mapper = Mapper(sequence)

    # state_vector = {}
    regionList = list(regionMap.keys())
    for state in stateList:

      for k in range(len(regionList)):

        region_name = regionList[k]
        region_prev = None if k == 0 else regionList[(k-1)]
        region_next = None if k == len(regionList)-1 else regionList[(k+1)]
        region = regionMap[region_name]

        stateRegion = {
          'score': 0,
          'start': -1,
          'stop': -1,
          'size': 0,
          'vector': []
        }

        # Add any eventual missing start position
        state_start = state[0] if len(state) > 0 else None
        if state_start and state_start['hmm'] > region['start']+1:
          for i in range(region['start']+1, min(region['stop']+1, state_start['hmm'])):
            stateRegion['vector'].append({'hmm': i, 'type': '-', 'idx': None, 'score': 0})
            # print('Adding {0} to {1}'.format(i, region_name))


        # Assemble region and calculate score
        stateSize = 0
        for i in range(len(state)):
          hmm_pos = state[i]['hmm']
          if hmm_pos <= region['start']:
            continue

          if hmm_pos > region['stop']:
            break

          stateRegion['size'] += 1 if state[i]['type'] == 'm' else 0
          stateRegion['score'] += state[i]['score']
          stateRegion['start'] = stateRegion['start'] if stateRegion['start'] and stateRegion['start'] >= 0 else state[i]['idx']
          stateRegion['stop'] = state[i]['idx'] if state[i]['idx'] else stateRegion['stop']
          stateRegion['vector'].append(state[i])

        # Add any eventual missing end position
        state_stop = stateRegion['vector'][-1] if len(stateRegion['vector']) > 0 else None
        if state_stop and (state_stop['hmm']+1) < region['stop']:
          for i in range(state_stop['hmm']+1, region['stop']+1):
            stateRegion['vector'].append({'hmm': i, 'type': '-', 'idx': None, 'score': 0})

        # # Check if the region is complete
        if stateRegion['size'] == 0:
          continue

        # Append to mapper
        mapper.append(region_name, stateRegion)
        
    return mapper.process()

  def __hmm_annotate__(self, state, sequence):

    annotationMap = {
      1: 'FR1',
      27: 'CDR1',
      39: 'FR2',
      56: 'CDR2',
      66: 'FR3',
      105: 'CDR3',
      118: 'FR4'
    }

    # Produce annotations
    annotationList = []
    current = {'type': None, 'start': 0, 'stop': 0}
    numbered_size = len(state)
    alignedSequence = ''
    for j in range(numbered_size): # Iterate over domains

      position = state[j][0][0]
      status = state[j][0][1]
      residue = state[j][1]
      alignedSequence += sequence[residue] if residue != None else '-'

      annotation_key = annotationMap.pop(position, None)
      if annotation_key:
        if current['type'] != None and current['start'] != None and current['stop'] != None:
          annotationList.append(current)

        current = {'type': annotation_key, 'start': residue, 'stop': residue}

      # Update current
      current['start'] = current['start'] if current['start'] else residue
      current['stop'] = residue if residue else current['stop']

    # Append last fragment
    if current['type'] != None and current['start'] != None and current['stop'] != None:
      annotationList.append(current)

    return alignedSequence, annotationList

  def __hmm_validate__(self, state, sequence):

    last = -1
    nseq=""
    for (_, status), index in state:

      if index == None:
        continue

      assert index >= last, "Numbering was found to decrease along the sequence %s. Please report."% sequence.getName()
      last = index
      nseq += sequence[index].replace('-', '')

    # if nseq not in seq.replace("-",""):
    if nseq not in sequence.getSequence().replace('-', ''):
      print(nseq)
      print(sequence.getSequence().replace('-', ''))
      raise Exception("The algorithm did not number a contiguous segment for sequence %s. Please report" % sequence.getName())

    return True

  @property
  def details(self):
    '''
    Returns the model details
    '''
    return self._dataset

  def size(self, species, chain):
    '''
    Returns the model size for a specific pair of species and chain names
    '''

    name = '{species}_{chain}'.format(species=species, chain=chain)
    if name not in self._dataset:
      raise Exception('Unable to find {0}-{1} in HMMR model'.format(species, chain))

    return self._dataset[name]['size']
    