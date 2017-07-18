# PyParanoid

**Ryan Melnyk, Haney Lab, UBC**  
**v1.0.0 - July 2017**
***ryan.melnyk@msl.ubc.ca***


## Installation

PyParanoid is primarily written in python.  To manage python installations and environments, I highly recommend using Miniconda (https://conda.io/miniconda.html)

The current version of PyParanoid (v1.0.0) depends on a modified version of the InParanoid script (Sonnhammer et al., 2015 - http://inparanoid.sbc.su.se/cgi-bin/faq.cgi) and thus also requires perl to be installed.

### Python modules
```
conda install biopython
conda install pandas
conda install seaborn
```

### Other Software

#### Using Homebrew

OSX users can install these dependencies easily using [Homebrew](https://brew.sh/) and running the following commands.

```
brew tap homebrew/science
brew install diamond
brew install hmmer
brew install mcl
```

#### Manual installation

If Homebrew doesn't work or you aren't in OSX, you will have to install the dependencies manually. Please ensure that all executables are located in a folder accessible in your $PATH.

##### DIAMOND
See http://github.com/bbuchfink/diamond for details.  Linux users can download an executable [here](https://github.com/bbuchfink/diamond/releases).

##### HMMER
http://hmmer.org/download.html - Download the newest version for your operating system.  PyParanoid should work with HMMER2 and HMMER3.

##### mcl
mcl can be found [here](https://www.micans.org/mcl/index.html?sec_software).

## Running PyParanoid

PyParanoid is a pipeline for rapid identification of homologous gene families in a set of reference genomes. These gene families are assembled into a database which is then used to annotate additional strains and inform homology-driven comparative genomics.




## cookbook

```
#basic functionality
PyParanoid.py
ClusterHomologs.py
ExtractClusters.py
```

```
#use profile HMMs to annotate many more genomes
PropagateGroups.py
MakeLocusTagMatrix.py
```

```
#QC of clustering
RarefactionPangenome.py
ClassifyAllSequences.py
PlotGroupSizes.py
```
