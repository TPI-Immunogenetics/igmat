# IgMAT: Antibody Multispecies Annotation Tool

IgMAT is a tool for the automatic discrimination and annotation of antibody sequences, specifically designed to be integrated into analysis pipelines or being used as a cli tool. IgMAT is highly customizable, allowing the addition of custom antibody sequences datasets and generating a range of output formats including a bed file of FR and CDR coordinates allowing simple downstream analysis on individual regions.

## Requirements ##
* Python version 3.8 or greater
* [Hmmer](http://hmmer.org/download.html)

## Installing ##
IgMAT can be installed locally or in a python environment: 

    git clone git@github.com:TPI-Immunogenetics/igmat.git
    cd igmat
    python3 -m venv env
    source env/bin/activate
    pip install ./

The default HMM model generated from IMGT data needs to be build in order to run IgMAT:

    igmat build

Once done, to exit the environment, type:

    deactivate
## Configuration ##
IgMAT will automatically generate a configuration file in the user home directory: 

    ~/.igmat/config.yaml

All generated hmmer models will be stored in this folder and the configuration file will be automatically initialised with the path containing the hmmer executables. For each one of the cli tools, the option `--hmmerpath` is available to temporarily override the hmmer path.

## CLI tools ##
IgMAT comes with a set of cli tools for handling custom HMM models and processing data:

    igmat <tool> --help

Where `tool` is one of the following:
 
 - **run**: process the input sequence/file and run the annotation script.
 - **build**: build custom HMM models 
 - **list**: handles the available HMM models

## Examples ##
IgMAT can be used as a stand alone tool, or embedded into a custom script. Please check the [examples](/docs/examples.md) and a tutorial about [embedding](/docs/embedding.md) IgMAT in your scripts.