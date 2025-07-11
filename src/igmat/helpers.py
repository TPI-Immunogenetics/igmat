import os
import sys
import platform
from pathlib import Path
from igmat import configs

def get_root():
  return Path(__file__).parent.parent

def get_dir_data():
  return os.path.join(os.path.expanduser('~'), '.igmat')

# TODO: move to package_data
# https://stackoverflow.com/questions/14211575/any-python-function-to-get-data-files-root-directory
def get_dir_binary(filename):

  # Try to guess the OS
  osname = platform.system().lower()
  if osname == 'linux' and 'microsoft' in platform.uname().release.lower():
    osname = 'win32'

  result = os.path.abspath(os.path.join(sys.prefix, 'muscle', osname, filename))
  print(result)
  if not os.path.exists(result):
    print('Unable to find binary for {filename} and OS {osname}. Fallback to local'.format(filename=filename, osname=osname))
    return filename
    # raise Exception('Unable to find binary for {filename} and OS {osname}'.format(filename=filename, osname=osname))

  return result

def get_alignment_format(input_path):
  '''
  This function check the format of the input alignment
  '''
  try:
    with open(input_path, 'r') as file:
      first_line = file.readline().strip()

      # Check for FASTA format
      if first_line.startswith('>'):
        return "fasta"

      # Check for Stockholm format
      if first_line == "# STOCKHOLM 1.0":
        return "stockholm"

    # If it matches neither format
    return None

  except FileNotFoundError:
    print(f"File {input_path} not found.")
    return None