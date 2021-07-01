import logging

def install(verbosity=0):
  logging.setLoggerClass(VerboseLogger)

class VerboseLogger(logging.Logger):

  levels = {
    0: logging.ERROR,
    1: logging.WARN,
    2: logging.INFO,
    3: logging.DEBUG,
  }

  def __init__(self, *args, **kw):
    logging.Logger.__init__(self, *args, **kw)
    self.parent = logging.getLogger()

  def setVerbosity(self, level=0):
    level = min(3, max(0,level))
    self.parent.setLevel(self.levels[level])
