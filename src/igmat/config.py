import os
import yaml

class Config:

  def __init__(self, path, default=dict()):
    self._path = path
    self._data = self.__load(default)

  def __load(self, default):

    # Check for config file
    config = {}
    if os.path.exists(self._path):
      with open(self._path, 'r') as handle:
        config = yaml.full_load(handle) 

    # Merge default values
    config = self.__merge(config, default)

    # Store if not exits
    if not os.path.exists(self._path):
      with open(self._path, 'w') as handle:
        yaml.dump(config, handle)

    return config

  def __merge(self, src, dst):
    for key, value in src.items():
      if isinstance(value, dict):
        node = dst.setdefault(key, {})
        self.__merge(value, node)
      else:
        dst[key] = value

    return dst

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
