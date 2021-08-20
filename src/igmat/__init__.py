
import os
from igmat.config import Config

# Create igmat folder
root_path = os.path.join(os.path.expanduser('~'), '.igmat')
if not os.path.exists(root_path):
  os.mkdir(root_path)

# Load config data
configs = Config(os.path.join(root_path, 'config.yaml'))