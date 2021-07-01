
class Alphabet:

  # Li: Reduction of protein sequence complexity by residue grouping (DOI: 10.1093/protein/gzg044)
  available = {
    'full': 'ACDEFGHIKLMNPQRSTVWY',
    'li8':  'SYEEYGNIKIINPEKSSIYY',
    'li10': 'SYEEYGNVKLLNPEKSSVYY',
    'li12': 'ACDEYGNVKLLNPEKTTVYY',
    'li14': 'ACDEYGHVKLLNPEKSTVYY',
    'li16': 'ACDEYGHVKLLNPERSTVWY',
    'li18': 'ACDEYGHVKLMNPQRSTVWY'
  }

  @staticmethod
  def list():
    return Alphabet.available.keys()

  @staticmethod
  def reduce(sequence, alphabet='full'):

    if alphabet == 'full':
      return sequence

    if alphabet not in available:
      raise Exception('Invalid alphabet: \'%s\'' % alphabet)

    result = ''
    special = ['.', '-', '*']
    for char in sequence:

      if char in special:
        result += char
        continue

      index = available['full'].find(char)
      if index < 0:
        raise Exception('Invalid character: \'%s\'' % char)

      # Append character
      result += available[alphabet][index]

    return result