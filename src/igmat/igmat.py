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

  return result