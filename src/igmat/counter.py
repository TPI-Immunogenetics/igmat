from multiprocessing import RawValue, Lock

class Counter(object):

  def __init__(self, valList):

    self._list = dict()
    self._lock = Lock()
    for valType in valList:
      self._list[valType] = RawValue('i', 0)

  def increment(self, valType):
    with self._lock:

      if valType not in self._list:
        raise Exception('Invalid value name: %s' % valType)

      self._list[valType].value += 1

  def value(self, valType):
    with self._lock:

      if valType not in self._list:
        raise Exception('Invalid value name: %s' % valType)

    return self._list[valType].value