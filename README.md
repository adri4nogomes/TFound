# TFound: A software to decipher TF-cis regulatory element interaction
TFound provides a simple visualization of possible co-locations between TFs, generating graphs and tables that facilitate the display and understanding of such information. TFound was developed in Python and its graphical interface developed implementing Tkinter and Matplotlib libraries. As proof of concept, we investigated the Saccharomyces cerevisiae regulatory machinery, loading its reference genome, annotation and a set of position matrices from The Yeast Transcription Factor Specificity Compendium (YeTFaSCo) collection. TFound is a modular platform and supports other annotated genomes and sets of Position Matrices, providing a simple and accessible tool to decipher multiple regulatory mechanisms through TF-cis regulatory element interaction inference, paving the way to understanding and engineering of a range of microorganisms.

## Instalation
<ul>
  <li>Install Conda (https://conda.io/docs/user-guide/install/index.html)</li>
  <li>Install dependencies:</li>
  <ul>
    <li>BeautifulSoup4 (conda install beautifulsoup4);</li>
    <li>Numpy (conda install numpy);</li>
    <li>Pandas (conda install pandas);</li>
    <li>MatPlotLib (conda install matplotlib);</li>
    <li>Numba (conda install numba).</li>
  </ul>
  <li>Run DbCreator.py to download the PWMs from YeTFaSCo and TF information from yeastgenome.org to create the database.</li>
  <li>Run InterfaceView.py</li>
</ul>
