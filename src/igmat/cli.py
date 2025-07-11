import os
import sys
import shutil
import argparse
import traceback

from .alphabet import Alphabet
from .build import build
from .run import run
from . import igmat
from . import helpers
from . import configs

# Pretty table
from prettytable import PrettyTable

from igmat.hmm.manager import Manager

def __run(args):

  run(
    input=args.input, 
    model=args.model, 
    restrict=args.restrict, 
    logPath=args.log, 
    annotationPath=args.annotation, 
    ncpu=args.ncpu, 
    bit_score_threshold=args.bit_score_threshold,
    verbose=args.verbose, 
    hmmerpath=args.hmmerpath
  )

def __build(args):

  # Generate the output path
  os.makedirs(helpers.get_dir_data(), exist_ok=True)

  # Execute the build script
  build(args.name, 
    input=args.input, 
    alignment_path=args.alignment, 
    alphabet=args.alphabet, 
    hmmerpath=args.hmmerpath, 
    verbose=args.verbose
  )

def __list(args):

  if args.clear and args.details:
    raise Exception('Please specify one single flag')

  manager = Manager(args.hmmerpath)

  # Remove data
  if args.clear:
    print('Clearing data')
    manager.clear(args.name)
    return

  # Show details
  if args.details:

    if not args.name:
      raise Exception('No model name is provided')

    print('Showing details for model {name}'.format(name=args.name))
    table = PrettyTable(['Species', 'chain', 'size', 'entropy'])
    details = manager.details(args.name)
    for key in details:
      table.add_row([
        details[key]['species'],
        details[key]['chain'],
        details[key]['size'],
        details[key]['relent']
      ])

    # Print the table
    print(table)
    return

  # No flags defined. show a list of libraries
  table = PrettyTable(['name', 'alphabet', 'chains'])
  for library in manager.list():

    table.add_row([
      library['name'],
      library['alphabet'],
      ','.join(library['chain'])
    ])

  print(table)

def main(*args, **kwargs):

  parser = argparse.ArgumentParser(prog="IgMAT", formatter_class=argparse.RawDescriptionHelpFormatter)
  subparsers = parser.add_subparsers(help='sub-command help', dest='cmd')
  
  # Run command
  parser_run = subparsers.add_parser('run', help='Run IgMAT')
  parser_run.add_argument('--sequence','-i', type=str, help="A sequence or an input fasta file", dest="input")
  parser_run.add_argument('--model', '-m', type=str, help='Model name [default=IMGT]', default='IMGT')
  parser_run.add_argument('--restrict','-r', type=str, nargs="+", choices=["ig","tr","heavy", "light", "H", "K", "L", "A", "B"], default=[], help="Only recognise certain types of receptor chains.", dest="restrict")    
  parser_run.add_argument('--annotation','-an', type=str, default=False, help="Output annotation file name", dest="annotation")
  parser_run.add_argument('--log','-l', type=str, default=False, help="Log annotation results", dest="log")
  parser_run.add_argument('--ncpu','-p', type=int, default=1, help="Number of parallel processes to use. Default is 1.", dest="ncpu")
  parser_run.add_argument('--bit_score_threshold', type=int, default=80, help="Change the bit score threshold used to confirm an alignment should be used.", dest="bit_score_threshold")
  parser_run.add_argument('--hmmerpath','-hp', type=str, default="", help="The path to the directory containing hmmer programs. (including hmmscan)", dest="hmmerpath")
  parser_run.add_argument('--verbose', '-v', action="store_true", dest="verbose", help="run in verbose mode")
  parser_run.set_defaults(cmd=__run)
  
  # Build command
  parser_build = subparsers.add_parser('build', help='Build IgMAT libraries')
  parser_build.add_argument('--input','-i', type=str, help="Input sequence folder", dest="input", required=False, default=None)
  parser_build.add_argument('--name','-n', type=str, help="HMM model name", dest="name", required=False, default='IMGT')
  parser_build.add_argument('--alignment', '-a', type=str, help="The alignment file in stockholm format to be used to generate the HMM model", default=None)
  parser_build.add_argument('--alphabet', '-l', type=str, help='Sequence alphabet', dest="alphabet", default='full', choices=Alphabet.list())
  parser_build.add_argument('--hmmerpath','-hp', type=str, default="", help="The path to the directory containing hmmer programs. (including hmmscan)", dest="hmmerpath")
  parser_build.add_argument('--verbose', '-v', action="store_true", dest="verbose", help="run in verbose mode")
  parser_build.set_defaults(cmd=__build)

  # List command
  parser_list = subparsers.add_parser('list', help='List IgMAT libraries')
  parser_list.add_argument('--name', '-n', type=str, help="HMM model name", dest="name", required=False, default=None)
  parser_list.add_argument('--clear', '-c', action="store_true", dest="clear", help="Remove all libraries")
  parser_list.add_argument('--details','-d', action="store_true", help='Show HMM model details', dest='details')
  parser_list.add_argument('--hmmerpath', '-hp', type=str, default="", help="The path to the directory containing hmmer programs. (including hmmscan)", dest="hmmerpath")
  parser_list.add_argument('--verbose', '-v', action="store_true", dest="verbose", help="run in verbose mode")
  parser_list.set_defaults(cmd=__list)

  # Check input arguments
  args = parser.parse_args()
  if not args.cmd:
    parser.print_help()
    sys.exit(1)

  try:

    # Check that hmmscan can be found in the path
    args.hmmerpath = args.hmmerpath if args.hmmerpath else configs.get('hmmerpath', '')
    hmmer_scan = os.path.join(args.hmmerpath, "hmmscan") if args.hmmerpath else shutil.which('hmmscan')
    if not args.hmmerpath and not (os.path.exists(hmmer_scan) and os.access(hmmer_scan, os.X_OK)):
      raise Exception("hmmscan was not found in the path. Either install and add to path or provide path with commandline option.")

    # Execute sub process
    args.cmd(args)
  except Exception as e:
    
    if args.verbose:
      traceback.print_exc()
      
    print(str(e))
    sys.exit(1)

if __name__ == "__main__":
  main()