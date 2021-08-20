import os
import sys
import platform
from pathlib import Path
from igmat import configs

def get_root():
  return Path(__file__).parent.parent

def get_dir_data():
  return os.path.join(os.path.expanduser('~'), '.igmat')
  # return os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')

# TODO: move to package_data
# https://stackoverflow.com/questions/14211575/any-python-function-to-get-data-files-root-directory
def get_dir_binary(filename):

  # Try to guess the OS
  osname = platform.system().lower()
  if osname == 'linux' and 'microsoft' in platform.uname().release.lower():
    osname = 'win32'

  # result = os.path.abspath(os.path.join(get_root(), 'bin', osname, filename))
  result = os.path.abspath(os.path.join(sys.prefix, 'muscle', osname, filename))
  if not os.path.exists(result):
    raise Exception('Unable to find binary for {filename} and OS {osname}'.format(filename=filename, osname=osname))

  return result