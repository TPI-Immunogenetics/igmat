import os
import sys
import shutil
import argparse

from .alphabet import Alphabet
from .build import build
from .run import run
from . import helpers


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
    alignment=args.alignment, 
    alphabet=args.alphabet, 
    hmmerpath=args.hmmerpath, 
    verbose=args.verbose
  )

def __list(args):

  print('Showing the list of parameters')

def main(*args, **kwargs):

  parser = argparse.ArgumentParser(prog="IgMAT", formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("-v", "--verbose",
                    action="store_true",
                    dest="verbose",
                    help="run in verbose mode")
  # parser.add_argument('--verbose', '-v', help='Set verbosity level', dest='verbose', action='count', default=0)
  parser.add_argument('--hmmerpath','-hp', type=str, default="", help="The path to the directory containing hmmer programs. (including hmmscan)", dest="hmmerpath")
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
  parser_run.set_defaults(cmd=__run)
  
  # Build command
  parser_build = subparsers.add_parser('build', help='Build IgMAT libraries')
  parser_build.add_argument('--input','-i', type=str, help="Input sequence folder", dest="input", required=False, default=None)
  parser_build.add_argument('--name','-n', type=str, help="HMM model name", dest="name", required=False, default='IMGT')
  parser_build.add_argument('--alignment', '-a', type=str, help="The alignment file in stockholm format to be used to generate the HMM model", default=None)
  parser_build.add_argument('--alphabet', '-l', type=str, help='Sequence alphabet', dest="alphabet", default='full', choices=Alphabet.list())
  parser_build.set_defaults(cmd=__build)

  # List command
  parser_list = subparsers.add_parser('list', help='List IgMAT libraries')
  parser_list.set_defaults(cmd=__list)

  # Check input arguments
  args = parser.parse_args()
  if not args.cmd:
    parser.print_help()
    sys.exit(1)

  # Check that hmmscan can be found in the path
  if args.hmmerpath:
    scan_path = os.path.join(args.hmmerpath, "hmmscan")
    if not (os.path.exists(scan_path) and os.access(scan_path, os.X_OK)):
      raise Exception("No hmmscan executable file found in directory: %s" % args.hmmerpath)
  elif not shutil.which("hmmscan"):
    raise Exception("hmmscan was not found in the path. Either install and add to path or provide path with commandline option.")

  # Execute sub process
  args.cmd(args)


if __name__ == "__main__":
  main()