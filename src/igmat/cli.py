import argparse

def main(*args, **kwargs):

  parser = argparse.ArgumentParser(prog="IgMAT", formatter_class=argparse.RawDescriptionHelpFormatter)
  subparsers = parser.add_subparsers(help='sub-command help')
  
  # Run command
  parser_run = subparsers.add_parser('run', help='Run IgMAT')
  parser_run.add_argument('--sequence','-i', type=str, help="A sequence or an input fasta file", dest="input")
  parser_run.add_argument('--model', '-m', type=str, help='Model name [default=IMGT]', default='IMGT')
  parser_run.add_argument('--verbose', '-v', help='Set verbosity level', dest='verbose', action='count', default=0)
  parser_run.add_argument('--restrict','-r', type=str, nargs="+", choices=["ig","tr","heavy", "light", "H", "K", "L", "A", "B"], default=False, help="Only recognise certain types of receptor chains.", dest="restrict")    
  parser_run.add_argument('--annotation','-an', type=str, default=False, help="Output annotation file name", dest="annotation")
  parser_run.add_argument('--log','-l', type=str, default=False, help="Log annotation results", dest="log")
  parser_run.add_argument('--hmmerpath','-hp', type=str, default="", help="The path to the directory containing hmmer programs. (including hmmscan)", dest="hmmerpath")
  parser_run.add_argument('--ncpu','-p', type=int, default=1, help="Number of parallel processes to use. Default is 1.", dest="ncpu")
  parser_run.add_argument('--bit_score_threshold', type=int, default=80, help="Change the bit score threshold used to confirm an alignment should be used.", dest="bit_score_threshold")
  parser_run.set_defaults(func=build)
  
  # Build command
  parser_build = subparsers.add_parser('build', help='Build IgMAT libraries')
  parser_build.set_defaults(func=build)

  # Check input arguments
  args = parser.parse_args()


if __name__ == "__main__":
  main()