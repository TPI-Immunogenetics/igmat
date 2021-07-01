
import os
import sys
import json
import tempfile
import traceback
import logging
from subprocess import Popen, PIPE

from . import fasta
from .alphabet import Alphabet
# from src.hmm.parser import Parser as hmmParser

# Import the HMMER parser from the distributed version of Biopython.
try: 
    from Bio.SearchIO.HmmerIO import Hmmer3TextParser as HMMERParser
    from Bio.SearchIO.HmmerIO import Hmmer3TextIndexer
except ImportError as e:
    print('Unable to import Biopython. Please install Biopython module.')
    sys.exit(1)

# This is the main class that handles HMM annotation
class HMMmodel():

  def __init__(self, path, model='IMGT', hmmerpath=None):

    self.modelname = model
    self.modelpath = os.path.join(path, model)
    self.hmmerpath = hmmerpath
    self.dataset = {}

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
      hmmstat = "hmmstat"
      if self.hmmerpath:
        hmmstat = os.path.join(self.hmmerpath, hmmstat)

      command = [ hmmstat, hmmerpath]
      process = Popen( command, stdout=PIPE, stderr=PIPE)

      # Parse the hmmstat output
      while True:
        line = process.stdout.readline().decode("utf-8")
        if not line:
          break

        if line[0] == '#':
          continue

        line = line.split()
        if len(line) < 5:
          continue

        name = line[1]
        size = int(line[5])
        self.dataset[name] = {
          'size': size
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

    # Create a fasta file for all the sequences. Label them with their sequence index
    # Todo: use tempfile.TemporaryDirectory
    fasta_filehandle, fasta_filename =  tempfile.mkstemp( ".fasta", text=True)
    with os.fdopen(fasta_filehandle,'w') as outfile:
      outfile.write(">%s\n%s\n" % (sequence.getName(), Alphabet.reduce(sequence.getSequence(), self.alpha) ))

    output_filehandle, output_filename = tempfile.mkstemp(".txt", text=True)

    # Run hmmer as a subprocess
    hmmscan = "hmmscan"
    if self.hmmerpath:
      hmmscan = os.path.join(self.hmmerpath, hmmscan)

    try:

      command = [ hmmscan, "--domT", "0", "-o", output_filename, self.modelpath + '.hmm', fasta_filename]
      process = Popen(command, stdout=PIPE, stderr=PIPE  )
      _, pr_stderr = process.communicate()
      if pr_stderr:
          raise Exception(pr_stderr)

      # parser = hmmParser()
      # parser.run(output_filename)
      # sys.exit(1)


      # Parse the result
      parser = Hmmer3TextIndexer(output_filename)
      query = parser.get(0)

      # Iterate over the matches of the domains in order of their e-value (most significant first)
      domainList = []
      hitList = []
      for hsp in sorted(query.hsps, key=lambda x: x.evalue):

        # Skipping matches different to the top one
        if hitList and hitList[0]['id'] != hsp.hit_id:
          continue

        # Check if in the exclude list
        species, ctype = hsp.hit_id.split('_')
        if restrict and ctype not in restrict:
          logging.info('skipping chain {chain}'.format(chain=ctype))
          continue

        # This is a valid domain
        domainList.append(hsp)

        # Store hit
        hitList.append({
          'id': hsp.hit_id,
          'description': hsp.hit_description,
          'evalue': hsp.evalue,
          'bitscore': hsp.bitscore,
          'bias': hsp.bias,
          'start': hsp.query_start,
          'end': hsp.query_end
        })

      # No valid match found
      if not domainList:
        return None

      # Generate a consensus of the found domains
      consensus = self.__hmm_align__(domainList, query.seq_len)

      # Validate the alignment
      self.__hmm_validate__(consensus['state'], sequence)

      consensus.update({
        'hits': hitList
        })

      return consensus
    finally:
      os.close(output_filehandle)

      # clean up
      os.remove(fasta_filename)
      os.remove(output_filename)
        
    # This shouldn't happen
    return None

  def germline(self, state_vector, sequence, chain_type, allowed_species=None):

    """
    Find the closest sequence identity match.
    """

    # Open the model germline file
    germlinePath = self.modelpath + '.json'
    if not os.path.exists(germlinePath):
      raise Exception('Unable to find germline file for model \'{model}\''.format(model=self.modelname))

    germline = {}
    with open(germlinePath, 'r') as handle:
      germline = json.load(handle)
    
    genes = {
      'v_gene': None,
      'j_gene': None
    }

    # Extract the positions that correspond to match (germline) states. 
    state_dict = dict( ((i, 'm'),None) for i in range(1,129))
    state_dict.update(dict(state_vector))
    state_sequence = ''.join([ sequence[state_dict[(i, 'm')]] if state_dict[(i,'m')] is not None else "-" for i in range(1,129) ])

    # Iterate over the v-germline sequences of the chain type of interest.
    # The maximum sequence identity is used to assign the germline 
    if chain_type in germline["V"]:
      if allowed_species is not None:
        if not all( [ sp in germline['V'][chain_type] for sp in allowed_species ] ): # Made non-fatal
          return None
      else:
        allowed_species = all_species

      resultList = []
      for species in allowed_species:

        if species not in germline["V"][ chain_type ]: 
          continue # Previously bug.

        for gene, germline_sequence in germline["V"][ chain_type ][ species ].items():
          resultList.append({
            'species': species,
            'gene': gene,
            'identity': get_identity(state_sequence, germline_sequence)
          })

      genes['v_gene'] = max(resultList, key=lambda x: x['identity'])
      
      # Use the assigned species for the v-gene for the j-gene. 
      # This assumption may affect exotically engineered abs but in general is fair.    
      species = genes['v_gene']['species']

      # TODO: cleanup
      if chain_type not in germline['J'] or species not in germline['J'][chain_type]:
        return genes

      resultList = []
      for gene, germline_sequence in germline["J"][ chain_type ][ species ].items():
        resultList.append({
          'species': species,
          'gene': gene,
          'identity': get_identity(state_sequence, germline_sequence)
          })
      genes['j_gene'] = max(resultList, key=lambda x: x['identity'])

    return genes

  def __hmm_align__(self, hspList, seq_length):

    regionMap = {
      'FR1': {'start': 0, 'stop': 26},
      'CDR1': {'start': 26, 'stop': 38},
      'FR2': {'start': 38, 'stop': 55},
      'CDR2': {'start': 55, 'stop': 65},
      'FR3': {'start': 65, 'stop': 104},
      'CDR3': {'start': 104, 'stop': 116},
      'FR4': {'start': 116, 'stop': 128}
    }

    stateList = []
    for index, hsp in enumerate(hspList):

      # Extract the strings for the reference states and the posterior probability strings     
      reference_string = hsp.aln_annotation["RF"]
      state_string = hsp.aln_annotation["PP"]

      # Extract the start an end points of the hmm states and the sequence
      # These are python indices i.e list[ start:end ] and therefore start will be one less than in the text file
      _hmm_start = hsp.hit_start
      _hmm_end = hsp.hit_end
       
      _seq_start = hsp.query_start
      _seq_end = hsp.query_end

      # Extact the full length of the HMM hit
      species, ctype = hsp.hit_id.split('_')
      _hmm_length = self.getModelSize(species, ctype)

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

      stateList.append(state)


    state_vector = {}
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
        state_stop = state[-1] if len(state) > 0 else None
        if state_stop and state_stop['hmm'] < region['stop']+1:
          for i in range(state_stop['hmm'], region['stop']+1):
            stateRegion['vector'].append({'hmm': i, 'type': '-', 'idx': None, 'score': 0})

        # Check if the region is complete
        if len(stateRegion['vector']) == 0 or stateRegion['vector'][0]['hmm']-1 > region['start'] or stateRegion['vector'][-1]['hmm'] < region['stop']:
          print('Incomplete region: {name}'.format(name=region_name))
          continue

        if stateRegion['size'] == 0:
          continue

        # Check if the region is continuous
        if region_name in state_vector:

          # Check if contiguous with prev region
          if (region_prev and region_prev in state_vector) and state_vector[region_prev]['stop'] > stateRegion['start']:
            # print('Non contiguous region with {region_prev}'.format(region_prev=region_prev))
            continue

          # Check if contiguous with next region
          if (region_next and region_next in state_vector) and stateRegion['stop'] > state_vector[region_next]['start']:
            # print('Non contiguous region with {region_next}'.format(region_next=region_next))
            continue


        # This is the best region so far
        if region_name not in state_vector or state_vector[region_name]['score'] < stateRegion['score']:
          prev_score = 0 if (region_name not in state_vector) else state_vector[region_name]['score']
          # print('Region: {region}; Previous score: {prev_score}; New score: {new_score}; start: {start}; stop: {stop}'.format(
          #   region=region_name, 
          #   prev_score=prev_score, 
          #   new_score=stateRegion['score'],
          #   start=stateRegion['start'],
          #   stop=stateRegion['stop'])
          # )

          state_vector[region_name] = stateRegion


    # state_vector = {}
    # regionList = list(regionMap.keys())
    # print(regionList)
    # # for region_name in regionMap:
    # for k in range(len(regionList)):

    #   region_name = regionList[k]
    #   region_prev = None if k == 0 else regionList[(k-1)]
    #   region_next = None if k == len(regionList)-1 else regionList[(k+1)]
    #   print(region_prev, region_name, region_next)

    #   # Align states to get the best regions
    #   region = regionMap[region_name]
    #   for state in stateList:

    #     stateRegion = {
    #       'score': 0,
    #       'start': -1,
    #       'stop': -1,
    #       'size': 0,
    #       'vector': []
    #     }

    #     # Add any eventual missing start position
    #     state_start = state[0] if len(state) > 0 else None
    #     if state_start and state_start['hmm'] > region['start']+1:
    #       for i in range(region['start']+1, min(region['stop']+1, state_start['hmm'])):
    #         stateRegion['vector'].append({'hmm': i, 'type': '-', 'idx': None, 'score': 0})

    #     # Assemble region and calculate score
    #     stateSize = 0
    #     for i in range(len(state)):
    #       hmm_pos = state[i]['hmm']
    #       if hmm_pos <= region['start']:
    #         continue

    #       if hmm_pos > region['stop']:
    #         break

    #       stateRegion['size'] += 1 if state[i]['type'] == 'm' else 0
    #       stateRegion['score'] += state[i]['score']
    #       stateRegion['start'] = stateRegion['start'] if stateRegion['start'] and stateRegion['start'] >= 0 else state[i]['idx']
    #       stateRegion['stop'] = state[i]['idx'] if state[i]['idx'] else stateRegion['stop']
    #       stateRegion['vector'].append(state[i])

    #     # Average the score by size
    #     # stateRegion['score'] = stateRegion['score']/max(1, stateRegion['size'])

    #     # Add any eventual missing end position
    #     state_stop = state[-1] if len(state) > 0 else None
    #     if state_stop and state_stop['hmm'] < region['stop']+1:
    #       for i in range(state_stop['hmm'], region['stop']+1):
    #         stateRegion['vector'].append({'hmm': i, 'type': '-', 'idx': None, 'score': 0})

    #     # Check if the region is complete
    #     if len(stateRegion['vector']) == 0 or stateRegion['vector'][0]['hmm']-1 > region['start'] or stateRegion['vector'][-1]['hmm'] < region['stop']:
    #       print('Incomplete region: {name}'.format(name=region_name))
    #       continue

    #     # Check if the region is continuous
    #     # if region_name in state_vector:

    #     #   # Check if contiguous with prev region

    #     #   if (region_prev and region_prev in state_vector):
    #     #     print('Previous region')
    #     #     print(state_vector[region_prev])

    #     #   if (region_next and region_next in state_vector):
    #     #     print('Next region')
    #     #     print(state_vector[region_next])

    #     #   if (region_prev and region_prev in state_vector) and state_vector[region_prev]['stop'] > stateRegion['start']:
    #     #     print('Non contiguous region with {region_prev}'.format(region_prev=region_prev))
    #     #     continue

    #     #   # Check if contiguous with next region
    #     #   if (region_next and region_next in state_vector) and stateRegion['stop'] > state_vector[region_next]['start']:
    #     #     print('Non contiguous region with {region_next}'.format(region_next=region_next))
    #     #     continue


    #     # This is the best region so far
    #     if region_name not in state_vector or state_vector[region_name]['score'] < stateRegion['score']:
    #       prev_score = 0 if (region_name not in state_vector) else state_vector[region_name]['score']
    #       print('Region: {region}; Previous score: {prev_score}; New score: {new_score}; start: {start}; stop: {stop}'.format(
    #         region=region_name, 
    #         prev_score=prev_score, 
    #         new_score=stateRegion['score'],
    #         start=stateRegion['start'],
    #         stop=stateRegion['stop'])
    #       )

    #       state_vector[region_name] = stateRegion

    # for i in state_vector:
    #   print(i, state_vector[i]['score'])
    #   for j in range(len(state_vector[i]['vector'])):
    #     print(state_vector[i]['vector'][j])

    # Check for non contiguous regions
    # This fixes problems when aligning two or more hmms and something is overlapping
    prev_valid = None
    prev_region = None
    for region_name in state_vector:
      
      first_valid = None
      last_valid = None
      for state in state_vector[region_name]['vector']:
        first_valid = first_valid if first_valid and first_valid['idx'] else state
        first_valid = first_valid if first_valid and first_valid['idx'] else None
        last_valid = state if state['idx'] else last_valid

      # This region is not contiguous, it needs to be removed
      if prev_valid and first_valid and (first_valid['idx']-prev_valid['idx'] > 1):
        
        region_remove = region_name if region_name in ['CDR1', 'CDR2', 'CDR3'] else prev_region
        if region_remove:
          state_vector[ region_remove ]['vector'] = []
          for i in range(regionMap[region_remove]['start'], regionMap[region_remove]['stop']):
            state_vector[ region_remove ]['vector'].append({
              'hmm': i+1,
              'type': '-',
              'idx': None,
              'score': 0
              })

      # Update for next
      prev_region = region_name
      prev_valid = last_valid if last_valid else prev_valid

    # Convert the result map to a tuple list
    # To be noted:
    # - hmm is one-based
    # - idx is zero-based
    result = []
    for key, value in state_vector.items():
      for state in value['vector']:
        result.append(((state['hmm'], state['type']), state['idx']))

    # Fill in the missing bits
    position=0
    state_start=-1
    state_end=0
    while position < len(result):

      # Update start and end position on the query sequence
      state_start = state_start if (state_start != None and state_start >= 0) else result[position][-1]
      state_end = result[position][-1] if result[position][-1] else state_end
      
      # Get position
      prev = { 'pos': position-1, 'hmm': result[position-1][0][0], 'idx': result[position-1][1], 'type': result[position-1][0][1]} if position > 0 else None
      curr = { 'pos': position, 'hmm': result[position][0][0], 'idx': result[position][1], 'type': result[position][0][1]}
      if curr['type'] != '-':
        position += 1
        continue

      # There is nothing before
      if prev is None or prev['idx'] is None:
        position += 1
        continue

      # This is an invalid position, try to find the next valid position
      last = None
      while position < len(result):
        last = { 'pos': position, 'hmm': result[position][0][0], 'idx': result[position][1], 'type': result[position][0][1]}
        if last['type'] != '-':
          break

        last = None
        position += 1

      if last is None:
        position += 1
        continue

      # print(last)
      # print(prev['idx']+1)
      if not last['idx']:
        continue

      idx = prev['idx']+1
      state_pos = curr['pos']
      hmm_pos = curr['hmm']
      while idx < last['idx']:

        if result[state_pos][0][1] != '-':
          result.insert(state_pos, ((hmm_pos, 'i'), idx))
          state_pos += 1
        else:
          result[state_pos] = ((hmm_pos, 'm'), idx)
          state_pos += 1
          hmm_pos += 1

        idx += 1

    return {
      'state': result,
      'start': state_start,
      'end': state_end
    }

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

  def getModelSize(self, species, chain):

    name = '{species}_{chain}'.format(species=species, chain=chain)
    if name not in self.dataset:
      raise Exception('Unable to find %s %s in HMMR model' % (species, chain))

    return self.dataset[name]['size']

# def get_identity(state_sequence, germline_sequence):
#   """
#   Get the partially matched sequence identity between two aligned sequences. 
#   Partial in the sense that gaps can be in the state_sequence.
#   """
#   # Ensure that the sequences are the expected length
#   assert len( state_sequence) == len(germline_sequence ) == 128
#   n, m = 0, 0
#   for i in range( 128 ):
#     if germline_sequence[i] == "-":
#       continue

#     if state_sequence[i].upper() == germline_sequence[i]: 
#       m+=1

#     n+=1

#   if not n:
#     return 0    

#   return float(m)/n
    
def annotate(
  sequence,
  model,
  threshold=0,
  restrict=[]):

  # Check if the sequence if valid
  if not sequence.isValid():
    raise AssertionError("Invalid sequence: %s" % sequence.getName())

  # Perform the alignments of the sequences to the hmm dataset
  result = model.run(sequence, threshold=threshold, restrict=restrict)
  if not result:
    raise Exception('Unable to find match')

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
  numbered_size = len(result['state'])
  alignedSequence = ''
  for j in range(numbered_size): # Iterate over domains

    position = result['state'][j][0][0]
    status = result['state'][j][0][1]
    residue = result['state'][j][1]
    if status != '-':
      alignedSequence += sequence[residue] if residue != None else '-'

    if (position in annotationMap):
      if current['type'] != None and current['start'] != None and current['stop'] != None:
        annotationList.append(current)

      current = {'type': annotationMap[position], 'start': residue, 'stop': residue}

    # Update current
    current['stop'] = residue if residue else current['stop']

  # Append last fragment
  if current['type'] != None and current['start'] != None and current['stop'] != None:
    annotationList.append(current)
  
  # Merge annotation results
  result.update({
    'alignment': alignedSequence,
    'annotation': annotationList,
  })

  return result