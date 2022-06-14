import os
import sys
import platform
import shutil
from pathlib import Path
from igmat import configs

def get_root():
  return Path(__file__).parent.parent

def get_os():

  # Try to guess the OS
  osname = platform.system().lower()
  if osname == 'linux' and 'microsoft' in platform.uname().release.lower():
    osname = 'win32'

  return osname

def get_dir_data():
  return os.path.join(os.path.expanduser('~'), '.igmat')
  # return os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')

def get_bin_path(name, path=None):
  '''
  Get the path of a binary file, adding the extensione if needed
  '''

  # Try to get the path from system
  bin_path = shutil.which(name)
  osname = platform.system().lower()
  if not bin_path and osname == 'win32':
    name += '.exe'

  if not os.path.exists(bin_path) and path is None:
    return None

  return bin_path if path is None else bin_path

  # Append extension for win32 systems
  # 
  # if osname == 'win32':
  #   name += '.exe'

  # return name if path is None else os.path.join(path, name)

# TODO: move to package_data
# https://stackoverflow.com/questions/14211575/any-python-function-to-get-data-files-root-directory
def get_dir_binary(filename):

  osname = get_os()
  result = os.path.abspath(os.path.join(sys.prefix, 'muscle', osname, filename))
  if not os.path.exists(result):
    print('Unable to find binary for {filename} and OS {osname}. Fallback to local'.format(filename=filename, osname=osname))
    return filename
    # raise Exception('Unable to find binary for {filename} and OS {osname}'.format(filename=filename, osname=osname))

  return result