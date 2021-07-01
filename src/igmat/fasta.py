import re

class sequence:
  def __init__(self, name, sequence):
    self.name = name
    self.sequence = sequence

  def __getitem__(self, arg):
    return self.sequence[arg]

  def getName(self):
    return self.name

  def getSequence(self):
    return self.sequence

  def getHash(self):
    return hash(self.sequence)

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