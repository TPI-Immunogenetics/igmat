import os
import json
import shutil
import igmat.helpers as helpers

from pathlib import Path
from igmat.hmm.model import Model

class Manager():

  def __init__(self, hmmerpath=None):
    self._hmmerpath = hmmerpath
    self._path = helpers.get_dir_data()

  def load(self, name):
    '''
    Load a model from data folder
    '''

    if not os.path.exists(os.path.join(self._path, name + '.json')):
      raise Exception('Unable to find model {name}'.format(name=name))

    return Model(self._path, name, self._hmmerpath)

  def list(self):
    '''
    Get a list of available models
    '''
    result = []
    if not os.path.exists(self._path):
      return result

    # Get the data directory
    for file in os.listdir(self._path):
      if not file.endswith('.json'):
        continue

      # Read the file
      with open(os.path.join(self._path, file)) as handle:
        data = json.load(handle)

      result.append({
        'name': data['name'],
        'alphabet': data['alphabet'],
        'chain': data['chain']
      })

    return result

  def clear(self, name=None):
    '''
    Remove all - or a specific - library
    '''

    def removeLibrary(path, name):

      # Removing library files
      shutil.rmtree(os.path.join(path, name), ignore_errors=True)
      for file in Path(path).glob("{name}.*".format(name=name)):
        file.unlink()

    # Name defined - remove a specific library
    if not os.path.exists(self._path):
      return

    if name:
      removeLibrary(self._path, name)
      return

    # Iterate the data directory
    for file in os.listdir(self._path):
      if not file.endswith('.json'):
        continue

      # Remove single library
      removeLibrary(self._path, file.replace('.json', ''))

  def details(self, name):
    '''
    Get details for a model
    '''

    model = self.load(name)
    # model = HMMmodel(helpers.get_dir_data(), name)

    return model.details