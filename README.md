# IgMAT: Antibody Multispecies Annotation Tool


Requirements:
 -  HMMER3 version 3.1b1 or higher - [http://hmmer.org/](http://hmmer.org/)
 - biopython [https://biopython.org/](https://biopython.org/)

e.g. 
    `ANARCI -i Example_sequence_files/12e8.fasta `
    This will number the files in 12e8.fasta with imgt numbering scheme and print to stdout.

    ANARCI -i Example_sequence_files/sequences.fasta -o Numbered_sequences.anarci -ht hit_tables.txt -s chothia -r ig 
    This will number the files in sequences.fasta with chothia numbering scheme only if they are an antibody chain (ignore TCRs).
    It will put the numbered sequences in Numbered_sequences.anarci and the alignment statistics in hit_tables.txt
    
    ANARCI -i Example_sequence_files/lysozyme.fasta
    No antigen receptors should be found. The program will just list the names of the sequences. 

    ANARCI -i EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSA
    Or just give a single sequence to be numbered. 

optional arguments:
  -h, --help            show this help message and exit
  --sequence INPUTSEQUENCE, -i INPUTSEQUENCE
                        A sequence or an input fasta file
  --model, -m <NAME>
                        The name of the HMMER model to use [default: IMGT]
  --verbose, -v
                        Set program verbosity
  --restrict {ig,tr,heavy,light,H,K,L,A,B} [{ig,tr,heavy,light,H,K,L,A,B} ...], -r {ig,tr,heavy,light,H,K,L,A,B} [{ig,tr,heavy,light,H,K,L,A,B} ...]
                        Restrict ANARCI to only recognise certain types of
                        receptor chains.
  --hmmerpath HMMERPATH, -hp HMMERPATH
                        The path to the directory containing hmmer programs.
                        (including hmmscan)
  --ncpu NCPU, -p NCPU  Number of parallel processes to use. Default is 1.
  --assign_germline     Assign the v and j germlines to the sequence. The most
                        sequence identical germline is assigned.
  --use_species {alpaca,rabbit,rhesus,pig,rat,human,mouse}
                        Use a specific species in the germline assignment.
  --bit_score_threshold BIT_SCORE_THRESHOLD
                        Change the bit score threshold used to confirm an
                        alignment should be used.

Author: James Dunbar (dunbar@stats.ox.ac.uk)
        Charlotte Deane (deane@stats.ox.ac.uk)

Contact: opig@stats.ox.ac.uk
        
--------------------------------------------------------------------------------------
# Setup #
In order to use Anarci, a HMMER dataset needs to be created. The utility script build.py generates a model starting either from the online IMGT reference dataset, or from a custom set of sequences.

- To build the IMGT HMMER model, please run `build.py` without any arguments. An updated list of sequences will be downloaded and a `.hmm` model will be generated.
 - To build a custom model, two sets of files need to be generated for each input species: a fasta file with sequences from the V-REGION, and a fasta file with sequences from the J-REGION. The files must be named in this format: `{species}_{chain}[V or J].fasta`, where *species* is the name of the organism, *chain* is the chain type (A,B,G,D,H,K,L). Additionally, a *name* parameter needs to be passed to the build script with the name of the model that will then be passed as input to the *anarci* script

Example of building a custom *cattle* dataset:
`build.py -i ./test/build/ -n cattle`

## Reduced alphabet ## 
The complexity in protein systems could be reduced by sorting natural occurring amino acids with similarities into groups, generating reduced alphabets [Li et al., 2003](https://doi.org/10.1093/protein/gzg044). This can be of help to improve performance with small datasets. To use this feature, the dataset need to be created with the `--alphabet <name>`, where `<name>` is the name of the alphabet (see below).
  
| Name | Number | g01 | g02 | g03 | g04 | g05 | g06 | g07 | g08 | g09 | g10 | g11 | g12 | g13 | g14 | g15 | g16 | g17 | g18 | g18 | g20
|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|
| li2 | 2 | CFYWMLIV | GPATSNHQEDRK 
| li3 | 3 | CFYWMLIV | GPATS | NHQEDRK 
| li4 | 4 | CFYW | MLIV | GPATS | NHQEDRK 
| li5 | 5 | CFYW | MLIV | G | PATS | NHQEDRK 
| li6 | 6 | CFYW | MLIV | G | P | ATS | NHQEDRK 
| li7 | 7 | CFYW | MLIV | G | P | ATS | NHQED | RK 
| li8 | 8 | CFYW | MLIV | G | P | ATS | NH | QED | RK 
| li9 | 9 | CFYW | ML | IV | G | P | ATS | NH | QED | RK 
| li10 | 10 | C | FYW | ML | IV | G | P | ATS | NH | QED | RK 
| li11 | 11 | C | FYW | ML | IV | G | P | A | TS | NH | QED | RK 
| li12 | 12 | C | FYW | ML | IV | G | P | A | TS | NH | QE | D | RK 
| li13 | 13 | C | FYW | ML | IV | G | P | A | T | S | NH | QE | D | RK 
| li14 | 14 | C | FYW | ML | IV | G | P | A | T | S | N | H | QE | D | RK 
| li15 | 15 | C | FYW | ML | IV | G | P | A | T | S | N | H | QE | D | R | K 
| li16 | 16 | C | FY | W | ML | IV | G | P | A | T | S | N | H | QE | D | R | K 
| li17 | 17 | C | FY | W | ML | IV | G | P | A | T | S | N | H | Q | E | D | R | K 
| li18 | 18 | C | FY | W | M | L | IV | G | P | A | T | S | N | H | Q | E | D | R | K 
| li19 | 19 | C | F | Y | W | M | L | IV | G | P | A | T | S | N | H | Q | E | D | R | K 
| full | 20 | C | F | Y | W | M | L | I | V | G | P | A | T | S | N | H | Q | E | D | R | K 

A dataset with the reduced alphabet will be created in the data folder, named `<name>_<alphabet>`, where `<name>` is the given model name, and `alphabet` is the alphabet name. By default, `full` is selected as no reduction in the alphabet. When `full` is selected, the alphabet name is omitted in the model name. 