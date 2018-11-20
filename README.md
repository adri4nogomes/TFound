# TFound: A software to decipher TF-cis regulatory element interaction
<p>TFound provides a simple visualization of possible co-locations between TFs, generating graphs and tables that facilitate the display and understanding of such information. TFound was developed in Python and its graphical interface developed implementing Tkinter and Matplotlib libraries. As proof of concept, we investigated the Saccharomyces cerevisiae regulatory machinery, loading its reference genome, annotation and a set of position matrices from The Yeast Transcription Factor Specificity Compendium (YeTFaSCo) collection. TFound is a modular platform and supports other annotated genomes and sets of Position Matrices, providing a simple and accessible tool to decipher multiple regulatory mechanisms through TF-cis regulatory element interaction inference, paving the way to understanding and engineering of a range of microorganisms.<p>

## Getting Started

### Prerequisites
<p>Download and install Conda (https://conda.io/docs/user-guide/install/index.html)</p>
<p>Install dependencies:</p>
<pre><code>conda install -y beautifulsoup4 numpy pandas matplotlib numba biopython</code></pre>

### Instalation
<p>Download the PWMs from YeTFaSCo and TF information from yeastgenome.org to create the database:</p>
<pre><code> python3 DbCreator.py</code></pre>

### Running
<pre><code>python3 InterfaceView.py</code></pre>
