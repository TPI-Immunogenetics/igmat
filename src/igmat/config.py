import os
import yaml

class Config:

  def __init__(self, path):
    self._path = path
    self._data = self.__load()

  def __load(self):

    # Check for config file
    config = {}
    if os.path.exists(self._path):
      with open(self._path, 'r') as handle:
        config = yaml.full_load(handle) 

    return config

  def get(self, name, default=None):

    """
    Use MongoDB style 'something.by.dot' syntax to retrieve objects from Python dicts.
    """
    val = self._data
    try:
      for key in name.split('.'):
        if key not in val:
          raise Exception('Invalid key')

        val = val[key]
    except:
      return default
    
    return val
