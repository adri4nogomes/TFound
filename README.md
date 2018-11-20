# TFound: A software to decipher TF-cis regulatory element interaction
<p>TFound provides a simple visualization of possible co-locations between TFs, generating graphs and tables that facilitate the display and understanding of such information. TFound was developed in Python and its graphical interface developed implementing Tkinter and Matplotlib libraries. As proof of concept, we investigated the Saccharomyces cerevisiae regulatory machinery, loading its reference genome, annotation and a set of position matrices from The Yeast Transcription Factor Specificity Compendium (YeTFaSCo) collection. TFound is a modular platform and supports other annotated genomes and sets of Position Matrices, providing a simple and accessible tool to decipher multiple regulatory mechanisms through TF-cis regulatory element interaction inference, paving the way to understanding and engineering of a range of microorganisms.<p>

<p align="center"><img src="https://github.com/adri4nogomes/TFound/blob/master/2.png" /></p>
<p align="center"><img src="https://github.com/adri4nogomes/TFound/blob/master/1.png"/></p>

## Getting Started

### Prerequisites
<ul>
  <li>Install Conda (https://conda.io/docs/user-guide/install/index.html)</li>
  <li>Run <code>conda install -y beautifulsoup4 numpy pandas matplotlib numba biopython</code> to install all dependencies.</li>
</ul>

### Instalation
<ul>
  <li>Clone this repo to your local machine.</li>
  <li><p>Run <code>python3 DbCreator.py</code> to download the PWMs from YeTFaSCo and TF information from yeastgenome.org to create the database.</p></li>
</ul>

## Running
<pre><code>python3 InterfaceView.py</code></pre>

## Author
<ul>
  <li><p><b>Adriano Gomes Silva </b> - <a href="http://silvarochar.wixsite.com/ssbl">Systems and Synthetic Biology Group </a></p></li>
</ul>

