# TFound: A software to decipher TF-cis regulatory element interaction
TFound provides a simple visualization of possible co-locations between TFs, generating graphs and tables that facilitate the display and understanding of such information. TFound was developed in Python and its graphical interface developed implementing Tkinter and Matplotlib libraries. As proof of concept, we investigated the Saccharomyces cerevisiae regulatory machinery, loading its reference genome, annotation and a set of position matrices from The Yeast Transcription Factor Specificity Compendium (YeTFaSCo) collection. TFound is a modular platform and supports other annotated genomes and sets of Position Matrices, providing a simple and accessible tool to decipher multiple regulatory mechanisms through TF-cis regulatory element interaction inference, paving the way to understanding and engineering of a range of microorganisms.

## Dependencies
<ul>
  <li>BeautifulSoup4 (https://www.crummy.com/software/BeautifulSoup/bs4/doc/);</li>
  <li>Numpy (http://www.numpy.org/ );</li>
  <li>Pandas (https://pandas.pydata.org/);</li>
  <li>MatPlotLib (https://matplotlib.org/);</li>
  <li>Numba (http://numba.pydata.org/).</li>
  <li>Biopython (https://biopython.org/)</li>
</ul>

## Instalation
<ul>
  <li>Install Conda (https://conda.io/docs/user-guide/install/index.html)</li>
  <li>Run: <i>conda install -y beautifulsoup4 numpy pandas matplotlib numba biopython</i>
  <li>Run: <i>python3 DbCreator.py</i> (to download the PWMs from YeTFaSCo and TF information from yeastgenome.org to create the database.)</li>
  <li>Run: <i>python3 InterfaceView.py</i></li>
</ul>
