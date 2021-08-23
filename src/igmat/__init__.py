
import os
import shutil
from igmat.config import Config

'''
Load the configuration file and return a dictionary
'''
def __init():

  # Create igmat folder
  root_path = os.path.join(os.path.expanduser('~'), '.igmat')
  if not os.path.exists(root_path):
    os.mkdir(root_path)

  # Extract the hmmer path from system
  hmmerpath = shutil.which('hmmscan')
  if hmmerpath: 
    hmmerpath = hmmerpath.replace('hmmscan', '')

  return Config(os.path.join(root_path, 'config.yaml'), {
    'hmmerpath': hmmerpath
  })

# Initialise the configuration
configs = __init()
