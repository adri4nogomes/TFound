# TFound: A software to decipher TF-cis regulatory element interaction
TFound provides a simple visualization of possible co-locations between TFs, generating graphs and tables that facilitate the display and understanding of such information. TFound was developed in Python and its graphical interface developed implementing Tkinter and Matplotlib libraries. As proof of concept, we investigated the Saccharomyces cerevisiae regulatory machinery, loading its reference genome, annotation and a set of position matrices from The Yeast Transcription Factor Specificity Compendium (YeTFaSCo) collection. TFound is a modular platform and supports other annotated genomes and sets of Position Matrices, providing a simple and accessible tool to decipher multiple regulatory mechanisms through TF-cis regulatory element interaction inference, paving the way to understanding and engineering of a range of microorganisms.

## Getting Started

### Prerequisites
Download and install Conda (https://conda.io/docs/user-guide/install/index.html)<br>
Download all dependencies:<br>
<code>conda install -y beautifulsoup4 numpy pandas matplotlib numba biopython</code>

### Instalation
Download the PWMs from YeTFaSCo and TF information from yeastgenome.org to create the database:<br>
<code> python3 DbCreator.py</code>

### Running
<code>python3 InterfaceView.py</code>
