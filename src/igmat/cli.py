import os
import sys
import argparse

from .alphabet import Alphabet
from .build import build
from . import helpers

def __build(args):

  file_path = helpers.get_binary('muscle')

  
  # Generate the output path
  os.makedirs(helpers.get_dir_data(), exist_ok=True)

  build(args.name, input=args.input, alignment=args.alignment, alphabet=args.alphabet, hmmerpath=args.hmmerpath, verbose=args.verbose)

def main(*args, **kwargs):

  parser = argparse.ArgumentParser(prog="IgMAT", formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('--verbose', '-v', help='Set verbosity level', dest='verbose', action='count', default=0)
  parser.add_argument('--hmmerpath','-hp', type=str, default="", help="The path to the directory containing hmmer programs. (including hmmscan)", dest="hmmerpath")
  subparsers = parser.add_subparsers(help='sub-command help')
  
  # Run command
  parser_run = subparsers.add_parser('run', help='Run IgMAT')
  parser_run.add_argument('--sequence','-i', type=str, help="A sequence or an input fasta file", dest="input")
  parser_run.add_argument('--model', '-m', type=str, help='Model name [default=IMGT]', default='IMGT')
  parser_run.add_argument('--restrict','-r', type=str, nargs="+", choices=["ig","tr","heavy", "light", "H", "K", "L", "A", "B"], default=False, help="Only recognise certain types of receptor chains.", dest="restrict")    
  parser_run.add_argument('--annotation','-an', type=str, default=False, help="Output annotation file name", dest="annotation")
  parser_run.add_argument('--log','-l', type=str, default=False, help="Log annotation results", dest="log")
  parser_run.add_argument('--ncpu','-p', type=int, default=1, help="Number of parallel processes to use. Default is 1.", dest="ncpu")
  parser_run.add_argument('--bit_score_threshold', type=int, default=80, help="Change the bit score threshold used to confirm an alignment should be used.", dest="bit_score_threshold")
  parser_run.set_defaults(cmd=__build)
  
  # Build command
  parser_build = subparsers.add_parser('build', help='Build IgMAT libraries')
  parser_build.add_argument('--input','-i', type=str, help="Input sequence folder", dest="input", required=False, default=None)
  parser_build.add_argument('--name','-n', type=str, help="HMM model name", dest="name", required=False, default='IMGT')
  parser_build.add_argument('--alignment', '-a', type=str, help="The alignment file in stockholm format to be used to generate the HMM model", default=None)
  parser_build.add_argument('--alphabet', '-l', type=str, help='Sequence alphabet', dest="alphabet", default='full', choices=Alphabet.list())
  parser_build.set_defaults(cmd=__build)

  # Check input arguments
  args = parser.parse_args()

  args.cmd(args)


if __name__ == "__main__":
  main()