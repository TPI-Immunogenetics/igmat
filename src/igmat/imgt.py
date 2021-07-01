
from html.parser import HTMLParser
from html.entities import name2codepoint

import os
import sys
import shutil
import traceback
import tempfile
import urllib.request
import urllib.parse
import urllib.error
from subprocess import Popen, PIPE

from . import fasta
from . import helpers


file_path  = os.path.split(__file__)[0]
bin_path = os.path.abspath(os.path.join(file_path, "../bin"))

class IMGTDBParser(HTMLParser):
  currenttag = None
  currentnamedent = None
  _data = []

  def handle_starttag(self, tag, attrs):
    self.currenttag=tag

  def handle_endtag(self, tag):
    self.currenttag=None

  def handle_data(self, data):
      split = data.split("\n")
      start = sum([ 1 if l[0]==">" else 0 for l in split if len(l)])
      if self.currenttag=="pre" and (self.currentnamedent ==">" or start):
          # Two different ways of parsing the html based on how IMGT have formatted the pages.
          # For some reason they format gene db differently sometimes (legacy?) 
          if start > 1: # If you encounter more than one line in the data with a fasta ">" symbol, all sequences will be in the same packet
              name, sequence = None, ""
              for l in split:
                  if not l: continue
                  if l[0]==">":
                      if sequence:
                          self._data.append( (name, sequence) )
                          name, sequence = None, ""
                      name = l
                  else:
                      sequence += l.replace(" ", "")
              if name and sequence:
                  self._data.append( (name, sequence) )
          else: # Otherwise it will be done entry by entry
              print("1")
              try:
                  name = split[0]
              except IndexError:
                  return
              sequence = ("".join( split[1:])).replace(" ", "")
              self._data.append( (name, sequence) )

  def handle_entityref(self, name):
      self.currentnamedent = chr(name2codepoint[name])

  def handle_charref(self, name):
      if name.startswith('x'):
          self.currentnamedent = chr(int(name[1:], 16))
      else:
          self.currentnamedent = chr(int(name))

  def rip_sequences(self,htmlstring):
      """
      Method for this subclass that automates the return of data
      """
      self.reset()
      self._data = []
      self.currenttag = None
      self.currentnamedent = None
      self.feed(htmlstring)
      self._data
      return self._data

def _mouse_delta(sequence):
  """
  Mouse delta chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.

  This is particularly bad because alignment is not even consistent within the chain and species!!!

  Remove and return
  """
  # Check in here because not all are bad...recheck again in the format v genes just to make sure.
  if sequence[103] != "C" or sequence[22] != "C":
      return sequence[ : 8 ] + sequence[ 9:85 ] + sequence[86:]
  return sequence

def _rhesus_lambda(sequence):
  """
  Rhesus lambda chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.
  Remove and return
  """
  return sequence[:20]+sequence[21:51]+ sequence[53:] 

def _mouse_alpha(sequence):
  """
  Mouse alpha chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.
  Remove and return
  """
  return sequence[:8]+sequence[9:85]+sequence[86:]

def _format_j_genes(jalignments):

  # The reference sequence 
  reference = {
    'name': 'Mus|H|Mus musculus|IGHJ3*01',
    'sequence': 'WFAYWGQGTLVTVSA',
    'start': 4,
    'stop': 19,
  }

  # Get alignment executable
  musclePath = helpers.get_dir_binary('muscle')
  # musclePath = os.path.join(bin_path, 'win32', 'muscle')
  # if sys.platform == 'linux' or sys.platform == 'linux2':
  #   musclePath = os.path.join(bin_path, 'linux', 'muscle')
  # elif sys.platform == 'darwin':
  #   musclePath = os.path.join(bin_path, 'osx', 'muscle')

  # musclePath = os.path.join(bin_path, "muscle")
  # if sys.platform == "darwin":
  #   musclePath = os.path.join(bin_path, "muscle_macOS")

  # Generate a temporary folder
  #with tempfile.TemporaryDirectory() as tmp:
  #  path = os.path.join(tmp, 'something')

  # Process sequences
  results = {}
  tempFolder = tempfile.mkdtemp()
  try:

    # Write a fasta file containing all sequences
    hasReference = False
    input_file = os.path.join(tempFolder, 'input.fa')
    output_file = os.path.join(tempFolder, 'output.fa')
    with open(input_file, "w" ) as outfile:
      for al in jalignments:
        for s in jalignments[al]:
          seqname = ">%s|%s|%s|%s" % tuple( list(al)+list(s) )
          outfile.write('>%s\n' % seqname)
          outfile.write('%s\n' % jalignments[al][s])

          if seqname == reference['name']:
            hasReference = True

      # If the reference sequence is not in the list, add it
      if not hasReference:
        outfile.write('>%s\n' % reference['name'])
        outfile.write('%s\n' % reference['sequence'])

    # Run muscle
    process = Popen( [ musclePath, "-in", input_file, "-gapopen", "-10", "-out", output_file, ], stdout=PIPE, stderr=PIPE )
    _, pr_stderr = process.communicate()
    #if pr_stderr:
    # raise Exception(pr_stderr)

    if not os.path.exists(output_file):
      raise Exception('Error while executing Muscle.')

    # Extract reference sequence
    for sequence in fasta.parse(output_file):
      if sequence.getName() == reference['name']:
        ref_aligned = sequence.getSequence()
        break

    start = ref_aligned.index(reference['sequence'])
    START = (start+1-reference['start']) if start > reference['start'] else 0
    END = start + len(reference['sequence'])
    for sequence in fasta.parse(output_file):

      if not hasReference and sequence.getName() == reference['name']:
        continue

      species, chain_type, id1, id2 = sequence.getName().strip(">").split("|")
      if (species, chain_type) not in results:   
          results[(species, chain_type)] = {}

      # We take the last 13 of the new alignment and pad into 20 long string 
      results[(species, chain_type)][ (id1, id2) ] = sequence.getSequence()[START: END][-14:].rjust(20).replace(" ", ".")

  #except Exception as e:
  #  print('Unable to format J genes: %s' % e)
  #  traceback.print_exc()
  finally:

    # Remove temp folder
    shutil.rmtree(tempFolder)

  return results

def _format_v_genes(valignments):

  # These are special functions for specific chains
  speciesVFormat = {
    'Macaca+mulatta_LV': _rhesus_lambda,
    #'rhesus_LV': _rhesus_lambda,
    'Mus_AV': _mouse_alpha,
    #'mouse_AV': _mouse_alpha,
    #'mouse_DV': _mouse_delta
    'Mus_DV': _mouse_delta
  }
  
  results = {} 
  for entry in valignments:

    species, chain_type = entry 
    results[entry] = {}
    for seq in valignments[entry]:

      try:

        sequence = valignments[entry][seq]
        #if species in translations:
        #name = '{species}_{type}{chain}'.format(species=translations[species], type=chain_type, chain='V')
        name = '{species}_{type}{chain}'.format(species=species, type=chain_type, chain='V')
        if name in speciesVFormat:
          #print('Applying specific formatting for %s %s chain' % (species, chain_type))
          sequence = speciesVFormat[name](sequence)

        # Trim sequence to the right size (108) and pad with gaps on the right side
        results[entry][seq] = sequence[:108].ljust(108).replace(" ",".")

      except Exception as e:

        print('Unable to format V genes: %s' % e)

  return results

def format(alignments):

  # Format V region
  if 'V' not in alignments:
    raise Exception('Unable to find V-regions in alignments')

  alignments['V'] = _format_v_genes(alignments['V'])

  # Format J region
  if 'J' not in alignments:
    raise Exception('Unable to find J-regions in alignments')

  alignments['J'] = _format_j_genes(alignments['J'])

  return alignments

def extractIMGTFasta(url, filename, retry=3):

  # Check if the file already exists
  if os.path.isfile(filename):
    return True # The file already exists

  # Download data from IMGT
  count = retry
  data = None
  while count > 0:
    try:
      handle = urllib.request.urlopen(url)
      data = handle.read()

    except Exception as e:
      print('Unable to retrieve url %s [retry %d/%d]' % (url, count, retry))
      traceback.print_exc()
      count -= 1
      time.sleep(0.1)
      data = None
      continue

    # Data downloaded successfully
    break;

  # Something happens
  if data is None:
    raise Exception('Unable to retrieve url {url}'.format(url=url))

  # Parse extracted html data
  parser = IMGTDBParser()
  sequences = parser.rip_sequences(data.decode('utf-8'))
  if not sequences:
    raise Exception('No sequences available for URL \'%s\'' % url)
  
  # All good, store fasta file
  with open(filename, "w") as outfile:
    for name, sequence in sequences:
      outfile.write('>%s\n' % name)
      outfile.write('%s\n' % sequence)

  return True
