import os
import re
import sys
import json
import traceback
from subprocess import Popen, PIPE

# from src import fasta
import igmat.imgt as imgt
import igmat.igmat as igmat
import igmat.fasta as fasta
import igmat.helpers as helpers
from igmat.alphabet import Alphabet

from igmat.hmm.manager import Manager

urls = { 
  "HV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+IGHV&species={species}",
  "HJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+IGHJ&species={species}",
  "KV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+IGKV&species={species}",
  "KJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+IGKJ&species={species}",
  "LV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+IGLV&species={species}",
  "LJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+IGLJ&species={species}",
  "AV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+TRAV&species={species}",
  "AJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+TRAJ&species={species}",
  "BV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+TRBV&species={species}",
  "BJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+TRBJ&species={species}",
  "GV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+TRGV&species={species}",
  "GJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+TRGJ&species={species}",
  "DV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+TRDV&species={species}",
  "DJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+TRDJ&species={species}"
}

# Species as of 04-12-14
# Species as of 02-06-16 - alpaca added
# These are retrieved for all the antibodies
all_species = [
  "Homo+sapiens",
  "Mus",
  "Rattus+norvegicus",
  "Oryctolagus+cuniculus",
  "Macaca+mulatta",
  "Sus+scrofa",
  "Vicugna+pacos",
  # "Bos+taurus",
  "Ovis+aries"  # Added 08-2020
]

all_tr_species = [
  "Homo+sapiens",
  "Mus",
]

# file_path  = os.path.split(__file__)[0]

def output_stockholm(all_sequences, path):

  # Output a minimal stockholm alignment file for all sequences. 
  with open( path, "w") as outfile:
    for species, chain_type in all_sequences:
      sequences = all_sequences[(species, chain_type)]
      l = len(list(sequences.values())[0])
      assert all( [1 if l == len(sequences[s]) else 0 for s in sequences]), "Not all sequences in alignment are the same length"
      
      ID = '{species}_{chain}'.format(species=species, chain=chain_type)
      outfile.write("# STOCKHOLM 1.0\n")
      outfile.write("#=GF ID %s\n" % ID)
  
      pad_length = max(list(map(len, list(sequences.keys()))))+1
      for s in sequences:
        outfile.write("{name} {sequence}\n".format(name=s.replace(" ", "_").ljust(pad_length), sequence=sequences[s].replace(".","-")))

      outfile.write("{end} {match}\n".format(end="#=GC RF".ljust(pad_length), match="x"*len(sequences[s])))
      outfile.write("//\n")

  return path

def getIMGTData(path, verbose=False):

  imgt_fields =  [
    "accession_number",
    "allele",  
    "species",  
    "functionality",  
    "region",  
    "start_and_end_positions_IMGT/LIGM-DB_accession_number", 
    "number_of_nucleotides", 
    "codon_start", 
    "number_nucleotides_added_in_5'_compared_IMGT/LIGM-DB", 
    "number_nucleotides_added_in_3'_compared_IMGT/LIGM-DB", 
    "number_nucleotides_to_correct_sequencing_errors", 
    "number_of_amino_acids", 
    "number_of_characters", 
    "partial",  
    "reverse"
  ]

  # For all V and J gene types (H,K,L,A,B,G,D) parse IMGT database and extract fasta files
  sequenceData = {}
  for gene_type in urls:
    for species in all_species:

      # We don't want TCRs for all organisms
      chain_type = gene_type[0]
      region_type = gene_type[1]
      if chain_type in "ABGD" and species not in all_tr_species: 
        continue

      try:

        filename = os.path.join(path, "{species}_{gene}.fasta".format(species=species.replace("+", "_"), gene=gene_type))
        if not imgt.extractIMGTFasta(urls[gene_type].format(species=species), filename):
          raise Exception('Unable to extract %s %s from IMGT' % (species, gene_type))

        # All done
        print("Parsed and saved %s %s" % (species, gene_type))
        if region_type not in sequenceData:
          sequenceData[region_type] = {}

        records = {}
        for sequence in fasta.parse(filename):

          # Parse sequence name
          fields = [x.strip() for x in sequence.getName().split("|")]
          fields = dict(list(zip(imgt_fields, fields)))

          # Filter out sequences
          isPartial = True if fields['partial'] else False
          isReverse = True if fields['reverse'] else False
          isGermline = True if fields["allele"].split("*")[-1].strip() == "01" else False
          try:

            # Check sequence validity
            if not sequence.isValid(gap=True):
              raise Warning('Invalid sequence')

            # Skip sequence with no accession number
            if fields['accession_number'] == 'None':
              raise Warning('Unknown accession number')

            # print(isPartial, isReverse, isGermline)
            # if fields['allele'] == 'IGHV1-198*01':
            #   print(fields)
            #   sys.exit(1)
            #   
            
            # Skip non functional, incomplete or reverse sequences
            if fields["functionality"] != "F" or isPartial or isReverse:
              raise Warning('Non-functional [functionality=\'{0}\'], partial [partial=\'{1}\'] or reverse [reverse=\'{2}\'] sequence'.format(
                fields['functionality'], 
                fields['partial'],
                fields['reverse']
              ))


            # TODO: find a better way to exclude extremely similar sequences
            # if read_all:
            #   pass
            # elif fields["allele"].split("*")[-1].strip()!="01": 
            #   raise Warning('Skipping allele %s' % fields['allele'])

            # All good, store sequence
            records[ (fields["species"], fields["allele"] ) ] = sequence.getSequence()
          
          except Warning as e:

            if verbose:
              print(' Warning: sequence {0}: {1}'.format(fields['allele'], str(e)))

        # Store sequence data
        sequenceData[region_type][(species, chain_type)] = records

      except Exception as e:
        traceback.print_exc()
        print('Error while extracting IMGT data: %s' % e)

  print('Formatting IMGT data')
  sequenceData = imgt.format(sequenceData)

  return sequenceData

def getCustomData(inputPath, modelPath, imgtModel, verbose=False):

  reference = {
    'name': 'Mus|H|Mus musculus|IGHJ3*01',
    'sequence': 'WFAYWGQGTLVTVSA',
    'size': 15
  }

  sequenceData = {}
  for file in os.listdir(inputPath):

    # Skip unwanted files
    if not file.endswith((".fa",".fasta",".fas")):
      continue

    name = os.path.splitext(file)[0];
    match = re.search("^(\w+)\_(\w)([JV])$", name)
    if not match:
      continue

    data = {
      'species': match.group(1),
      'chain': match.group(2),
      'type': match.group(3)
    }

    print('Found chain %s of type %s from %s' % (data['chain'], data['type'], data['species']))

    # Generate a new entry for this species
    if data['type'] not in sequenceData:
      sequenceData[ data['type'] ] = {};

    # Parse sequence
    dataList = {}
    sequenceHash = []
    for sequence in fasta.parse(os.path.join(inputPath, file)):

      # Check for valid sequences
      if not sequence.isValid():
        raise Exception('Invalid characters in sequence \'{name}\''.format(name=sequence.getName()))
    
      if sequence.getHash() in sequenceHash:
        if verbose:
          print('Found duplicate sequence: {name}'.format(name=sequence.getName()))
        continue

      name = sequence.getName() \
        .replace('\t', '_') \
        .replace('*', '_') \
        .replace('|', '_')

      # Make sure to create an unique name
      count = 1
      seqName = name
      while (data['species'], seqName) in dataList:
        seqName = '{name}_{count}'.format(name=name, count=count)
        count += 1

      dataList[(data['species'], seqName)] = sequence.getSequence()
      sequenceHash.append(sequence.getHash())
      # dataList[(data['species'], name)] = sequence.getSequence()

    sequenceData[ data['type'] ][ (data['species'], data['chain'])] = dataList

  # Format V regions
  for chain in sequenceData['V']:

    for sequence in sequenceData['V'][chain]:

      # In order to get a nice result, we need to add a J-Region to the V-Region.
      # We are using a reference J-Region
      testSequence = sequenceData['V'][chain][sequence] + reference['sequence']
      result = igmat.annotate(fasta.sequence('Test', testSequence), imgtModel)
      # result = anarci.anarci_thread(fasta.sequence('Test', testSequence), imgtModel)
      if not result:
        continue

      # Now we need to remove the reference sequence, including any eventual gap
      referenceStop = result.start
      resultStop = 0 
      while referenceStop < len(sequenceData['V'][chain][sequence]):
        char = result.sequence[resultStop]
        if char == '-':
          resultStop += 1
          continue

        if char == sequenceData['V'][chain][sequence][referenceStop]:
          resultStop += 1
          referenceStop += 1

      # Trim out the reference sequence
      sequenceData['V'][chain][sequence] = result.sequence[:resultStop]

  # Apply IMGT formatting
  sequenceData = imgt.format(sequenceData)

  return sequenceData;

def getAlignmentData(path):

  sequenceData = {
    'V': {},
    'J': {}
  }
  try:
    with open(path, 'r') as handle:

      count = 1
      type = None
      for line in handle:

        # Get alignment type
        line = line.strip()
        if line.startswith('#=GF ID '):
          type = dict(list(zip(['species', 'chain'], line[8:].strip().split('_'))))

        # Skip comment line
        if not line or line[0] == '#':
          continue

        # End of alignment
        if line.startswith('//'):
          type = None
          continue

        # Extract sequence data
        line = line.split()
        if len(line) != 2:
          raise Exception('Invalid format at line %d' % count)

        head = line[0]
        sequence = line[1]
        if not type:
          raise Exception('Invalid or missing alignment ID. Must be in the format \'#=GF ID <species>_<chain>\'')

        # Extracting sequence name
        head = head.split('|')
        if len(head) != 3 or head[0] != type['species']:
          raise Exception('Invalid sequence name. Must be in the format \'<species>|<v-region>|<j-region>\'')

        # Validating sequence
        if len(sequence) != 128:
          raise Exception('Invalid sequence length. Must be exaclty 128 residues')

        if not (type['species'], type['chain']) in sequenceData['V']:
          sequenceData['V'][(type['species'], type['chain'])] = {}
          sequenceData['J'][(type['species'], type['chain'])] = {}

        # Add V region if not already existing
        if not (type['species'], head[1]) in sequenceData['V'][(type['species'], type['chain'])]:
          sequenceData['V'][(type['species'], type['chain'])][ (type['species'], head[1]) ] = sequence[:108]

        # Add J region if not already existing
        if not (type['species'], head[2]) in sequenceData['J'][(type['species'], type['chain'])]:
          sequenceData['J'][(type['species'], type['chain'])][ (type['species'], head[2]) ] = sequence[-20:]

        count += 1
  except Exception as e:
    print('  Error loading alignment at line %d: %s' % (count, e))
    return None

  return sequenceData

def generateAlignment(sequences, path, isIMGT=False, alph='full'):

  if not sequences:
    raise Exception('Empty alignment data')

  results = {}
  for species, chain_type in sequences['V']:
    
    if (species, chain_type) not in sequences['J'] or (species, chain_type) not in sequences['V']:
      continue

    # Do a pairwise combination of the v and j sequences to get putative germline sequences for the species.
    combined_sequences = {}
    for v in sequences['V'][ (species, chain_type) ]:
      vspecies, vallele = v

      # TODO: find a better way to exclude similar sequences
      if isIMGT and vallele.split("*")[-1].strip()!="01":
        print('Skipping sequence %s' % vallele)
        continue

      for j in sequences['J'][ (species, chain_type) ]:
        _, jallele = j

        # TODO: find a better way to exclude similar sequences
        if isIMGT and jallele.split("*")[-1].strip()!="01":
         continue

        # Combine the sequence
        combined = sequences['V'][ (species, chain_type) ][v] + sequences['J'][ (species, chain_type) ][j]

        # Translate to alphabet
        combined = Alphabet.reduce(combined, alph)

        # Store it
        name = ("%s|%s|%s"%(vspecies, vallele,jallele)).replace(" ", "_")
        combined_sequences[name] = combined

    results[ (species, chain_type) ] = combined_sequences

  # Write just the V and J combinations
  output_stockholm(results, path);

def validateAlignment(path, verbose=False):

  # ruleList = {
  #   23: ['C'],
  #   41: ['W'],
  #   89: ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V', 'P'],
  #   104: ['C'],
  #   118: ['F', 'W']
  # }
  # 

  ruleList = {
    23: {
      'residues': ['C'],
      'critic': True
    },
    41: {
      'residues': ['W'],
      'critic': False
    },
    89: {
      'residues': ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V', 'P'],
      'critic': False
    },
    104: {
      'residues': ['C'],
      'critic': True
    },
    118: {
      'residues': ['F', 'W'],
      'critic': False
    }
  }

  try:

    warningCount = 0
    if not os.path.exists(path):
      raise Exception('Unable to find alignment file')

    filename = os.path.basename(path)
    with open(path, 'r') as handle:

      count = 0
      for line in handle:

        # Skip comment line
        count += 1
        line = line.strip()
        if not line or line[0] == '#':
          continue

        # End of alignment
        if line.startswith('//'):
          type = None
          continue

        # Extract sequence data
        line = line.split()
        if len(line) != 2:
          raise Exception('Invalid format at line %d' % count)

        head = line[0]
        sequence = line[1]
        for index in ruleList:

          try:
                    
            if sequence[index-1] not in ruleList[index]['residues']:

              message = 'Line {0}: position {1} must be one of \'{2}\''.format(
                count,
                index,
                ','.join(ruleList[index]['residues'])
                )
              raise (Exception(message) if ruleList[index]['critic'] else Warning(message))
              # raise exception
              # if ruleList[index]['critic']:
              #   raise Exception('')
              # raise Warning('Position {0} must be one of: {1}'.format(
              #   index,
              #   ','.join(ruleList[index])
              #   ))

          except Warning as w:

            warningCount += 1
            if verbose:
              print('Warning {0}: {1}'.format(filename, w))
              print('  ' + sequence)
              print('  ' + (' ' * (index-1) + '*'))
            

  except Exception as e:
    print('Error validating alignment {0}: {1}'.format(filename, str(e)))
    return False

  if warningCount:
    print('Found {0} warnings while validating data. Please check the alignment'.format(warningCount))

  return True

def generateModel(alignment, name, path, hmmerpath=""):

  HMMPath = os.path.join(path, name + '.hmm')
  try:

    # Check for aligment path
    if not os.path.exists(alignment):
      raise Exception('Unable to find input alignment for hmmbuild')

    # Run hmmer as a subprocess
    hmmbuild = os.path.join(hmmerpath, "hmmbuild") if hmmerpath else "hmmbuild"
    command = [ hmmbuild, "--hand", HMMPath, alignment]
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    _, pr_stderr = process.communicate()

    if pr_stderr:
        raise Exception(pr_stderr)

    # Turn the output HMMs file into a binary form.
    hmmpress = os.path.join(hmmerpath, "hmmpress") if hmmerpath else "hmmpress"
    command = [ hmmpress, "-f", HMMPath]
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    _, pr_stderr = process.communicate()
    if pr_stderr:
       raise Exception(pr_stderr)

  except Exception as e:
    print('  Error building HMM: %s' % e)
    return False

  return True

def generateDatafile(name, sequences, filename, alph):

  data = {
    'name': name,
    'alphabet': alph,
    'chain': [],
    'J': {},
    'V': {}
  }

  # Store sequence data
  for chain in ['V', 'J']:
    for species, chain_type in sequences[chain]:
      for ((_,gene), seq) in sequences[chain][ (species, chain_type) ].items():

        # Add chain to output
        if chain_type not in data['chain']:
          data['chain'].append(chain_type)

        # Add chain key
        if chain_type not in data[chain]:
          data[chain][chain_type] = {}

        # Add species key
        if species not in data[chain][chain_type]:
          data[chain][chain_type][species] = {}

        # Add sequence
        if chain == 'V':
          assert len(seq)==108, species+_+gene+chain_type+_+seq+str(len(seq))
          data[chain][chain_type][species][gene] = seq.replace(".","-") + "-"*20
        elif chain == 'J':
          assert len(seq)==20
          data[chain][chain_type][species][gene] = "-"*108 + seq.replace(".","-")
        
  # Write to file
  with open(filename, 'w') as outfile:
    json.dump(data, outfile, indent=2)

def build(name, input=None, alignment=None, alphabet='full', hmmerpath=None, verbose=False):

  # Check input folder
  if input is not None:

    if not os.path.exists(input) or not os.path.isdir(input):
      raise Exception('The input argument must be a valid directory!')

    if alignment is not None:
      raise Exception('Please provide either an alignment file or an input folder')

  # Check model name
  if (input is not None or alignment is not None) and name == 'IMGT':
    raise Exception('Please provide a different name for a custom model')

  # Generate the output path
  output = helpers.get_dir_data()
  
  # Generate the model folder
  modelName = name + ('_' + alphabet if alphabet != 'full' else '')
  modelPath = os.path.join(output, modelName)
  os.makedirs(modelPath, exist_ok=True)

  # Create the manager object
  manager = Manager(hmmerpath)

  # Extract data from the provided input
  sequenceData = None
  isIMGT = True if name == 'IMGT' else False
  if isIMGT:
    sequenceData = getIMGTData(modelPath, verbose)
  else:

    # Generate the alignment from input data
    if input:

      # Load the IMGT HMMER model
      # try:
      # dataset = igmat.HMMmodel(output, 'IMGT', hmmerpath)
      dataset = manager.load('IMGT')
      # except Exception:
        # print('Unable to find IMGT HMM model. Please run the build script without arguments first.')
        # sys.exit()

      sequenceData = getCustomData(input, modelPath, dataset, verbose)

    # Load the alignment file
    elif alignment:

      sequenceData = getAlignmentData(alignment)

  # try:

  # Generate an alignment with all the combinations of V and J regions
  alignPath = os.path.join(modelPath, modelName + '.stockholm')
  print('Generating alignment data')
  generateAlignment(sequenceData, alignPath, isIMGT, alph=alphabet)

  # Validate the generated alignment
  if not validateAlignment(alignPath, verbose):
    raise Exception('Invalid alignment data')

  # Generate an HMM model with all the aligned sequences
  print('Generating HMM model \'%s\'' % modelName)
  if not generateModel(alignPath, modelName, output):
    raise Exception('Unable to generate HMM model')

  # Compile a file with the dataset properties and a dictionary containing all the v and j germline sequences.
  # germlinePath =  os.path.join(output, args.name + '.json')
  germlinePath = os.path.join(output, modelName + '.json')
  print('Generating data file \'%s\'' % modelName)
  generateDatafile(modelName, sequenceData, germlinePath, alph=alphabet)

  # except Exception as e:
  #   print('Error - %s' % e)
  #   sys.exit(1)