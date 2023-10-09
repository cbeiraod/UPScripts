# UniProt Tools
This repository started as an evolution of [SysBioTK](https://sysbiotk.sourceforge.io/). UniProt has updated and modernized the programmatic interface to their database since SysBioTK was created, so the need to adapt to this change was needed as well. The opportunity was taken in order to improve the library design and create a new set of tools.

The library is being reimplemented incrementally, with the functionality being implemented as needed.

The tools have been tested with python3 and are only guaranteed to work on this version of python.

Currently only the GO analysis tools are implemented and a script interface exists to run the tools on an individual basis.

**Scripts**:
  * go_ana.py - This script performs a GO analysis. The tool needs to be provided with the go.obo database file, which can be downloaded from https://geneontology.org/docs/download-ontology/. The tool should also be provided with a directory containing the XML output from the "ID Mapping" tool of the UniProt website saved as listUP.xml. The tool will create a summary of the GO terms associated with each UniProt Accession number and then will also produce a summary of how many proteins were tagged with each term, using a GO Slim (the default is Generic GO Slim) to reduce the amount of GO terms to be considered. If a protein is tagged with a GO term which is marked as being an "is a" or "part of" another GO term, this tree of relationships is parsed in order to find the filtered GO Slims each protein is tagged with.

## Usage
Please use a virtual environment (venv) to ensure all needed dependencies are installed and that there is no conflicts with the requirements from other scripts or the base system.

In a first step the venv needs to be created and the dependencies installed.
Then, every time before using UPTools, make sure that the venv is activated.

### Creation
Create the venv using the command: `python3 -m venv venv`

In the example above the venv is created in a directory called venv (the last parameter), feel free to change this. Remember to adapt subsequent commands to use your chosen venv name.

Active the venv with: `source venv/bin/activate`

Install the library dependencies with:
```
python -m pip install pillow
python -m pip install openpyxl
```

### Activation
Active the venv with: `source venv/bin/activate`

### Deactivation
Once you are done and if you no longer need the venv, it can be deactivated with: `deactivate`

## Developers
To be written... I want to implement this correctly as a library, perhaps even on pypi for ease of use and installation. Also use unit testing to ensure things are working as expected.
