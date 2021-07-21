# IgMAT: Antibody Multispecies Annotation Tool

IgMAT is a tool for the automatic discrimination and annotation of antibody sequences, specifically designed to be integrated into analysis pipelines or being used as a cli tool. IgMAT is highly customizable, allowing the addition of custom antibody sequences datasets and generating a range of output formats including a bed file of FR and CDR coordinates allowing simple downstream analysis on individual regions.

## Requirements ##
* Python version 3.8 or greater

## Installing ##
IgMAT can be installed locally or in a python environment: 

    python3 -m venv env
    source env/bin/activate
    pip install ./

The default HMM model generated from IMGT data needs to be build in order to run IgMAT:

    igmat build

Once done, to exit the environment, type:

    deactivate

## CLI tools ##
IgMAT comes with a set of cli tools for handling custom HMM models and processing data:

    igmat <tool> --help

Where `tool` is one of the following:
 
 - **run**: process the input sequence/file and run the annotation script.
 - **build**: build custom HMM models 
 - **list**: handles the available HMM models
