import os
from typing import Union

import igmat.fasta as fasta

from igmat.hmm.manager import Manager
from igmat.hmm.model import Model

def annotate(
  sequence: Union[str, fasta.sequence],
  model: Union[str, Model] = 'IMGT',
  threshold=0,
  restrict=[],
  hmmerpath=None
  ):

  # Generate a sequence object
  if isinstance(sequence, str):
    sequence = fasta.sequence('input', sequence)

  # Load the model object
  if isinstance(model, str):
    manager = Manager(hmmerpath)
    model = manager.load(model)

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