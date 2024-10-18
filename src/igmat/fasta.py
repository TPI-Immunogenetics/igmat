import re
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna

NUCLEOTIDE_CHARS = {'A', 'C', 'T', 'G', 'U', 'a', 't', 'c', 'g', 'u'}
PROTEIN_CHARS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

class sequence:
  
  def __init__(self, name, sequence):
    self.name = name
    self.sequence = sequence
    self.type = None

  def __getitem__(self, arg):
    return self.sequence[arg]

  def getName(self):
    return self.name

  def getSequence(self):
    return self.sequence

  def getSize(self):
    return len(self.sequence)

  def getHash(self):
    return hash(self.sequence)
  
  def getType(self):
    if self.type is not None:
      return self.type
    
    sequence = self.sequence.strip()
    try:
      # If all characters are valid nucleotide bases, classify as nucleotide
      if all(char in NUCLEOTIDE_CHARS for char in sequence):
        self.type = "nucleotide"

      # If all characters are valid protein residues, classify as protein
      elif all(char in PROTEIN_CHARS for char in sequence):
        self.type = "protein"
      else:
        raise Exception('Unknown protein type')
    except Exception as e:
      print(str(e))
      self.type = 'unknown'
      
    return self.type

  def isValid(self, gap=False):
    return True if re.search('[^QWERTYIPASDFGHKLCVNMX' + ('\.' if gap else '') + ']', self.sequence) is None else False

def parse(path):

  with open(path) as fp:
    name, seq = None, []
    for line in fp:
      line = line.rstrip()
      if line.startswith(">"):

        # Append previous
        if name: 
            yield sequence(name, ''.join(seq))

        name, seq = line[1:], []
      else:
        seq.append(line)

    # Append last in list
    if name: 
      yield sequence(name, ''.join(seq))
      
def translate(sequence, frameshift=None):
  
  allowed_frames = [-3, -2, -1, 1, 2, 3 ]
  if frameshift is None:
    frameshift = allowed_frames
    
  if not all(item in allowed_frames for item in frameshift):
    raise Exception('Invalid frameshift')
  
  result = {}
  
  # Define the forward and reverse complement sequences
  fwd_sequence = Seq(sequence)
  rvs_sequence = fwd_sequence.reverse_complement()
  for frame in frameshift:
    shift = abs(frame)-1
    if frame > 0:
      result[frame] = str(fwd_sequence[shift:].translate(to_stop=False))
    else:
      result[frame] = str(rvs_sequence[shift:].translate(to_stop=False))
      
  return result