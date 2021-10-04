
class Result():

  def __init__(self, sequence, start, end, annotation, hitList=None):
    self._sequence = sequence
    self._start = start
    self._end = end
    self._list = annotation
    self._hitList = hitList

  @property
  def sequence(self):
    return self._sequence

  @property
  def start(self):
    return self._start

  @property
  def end(self):
    return self._end

  def annotations(self):
    '''
    Returns a generator with the annotations
    '''

    for i in range(len(self._list)):
      yield {
        'type': self._list[i]['type'], 
        'start': self._list[i]['start'], 
        'stop': self._list[i]['stop'], 
      }

  def __str__(self):

    ''' 
    Generate a string representing the aligned sequence with separators 
    ''' 
    fragmentMap = {}
    for fragment in self._list:
      coordinates = (fragment['stop']-self._start)
      fragmentMap[coordinates] = fragment['type']

    sequence = ''
    count = 0
    header = ''
    for char in self._sequence:
      sequence += char
      if char == '-':
        continue

      if count in fragmentMap:
        
        header += fragmentMap[count].center(len(sequence)-len(header), '-') + '|'
        sequence += '|'

      count += 1

    # Generate the result header
    result = ''
    for hit in self._hitList:
      result += ' HMMR match: {} [eValue: {:.3e}]\n'.format(hit['id'], hit['evalue'])
          
    if result:
      result += '\n'

    # Append alignment data
    linesize = 70
    for i in range(0, len(sequence), linesize):
      result += ' {0:3d} {1} {2:3d}'.format(i+self._start, header[i:i+linesize], self._start+min(i+linesize, len(sequence))) + '\n'
      result += ' {0:3d} {1} {2:3d}\n'.format(i+self._start, sequence[i:i+linesize], self._start+min(i+linesize, len(sequence))) + '\n'

    return result