# TFound: A software to decipher TF-cis regulatory element interaction
TFound provides a simple visualization of possible co-locations between TFs, generating graphs and tables that facilitate the display and understanding of such information. TFound was developed in Python and its graphical interface developed implementing Tkinter and Matplotlib libraries. As proof of concept, we investigated the Saccharomyces cerevisiae regulatory machinery, loading its reference genome, annotation and a set of position matrices from The Yeast Transcription Factor Specificity Compendium (YeTFaSCo) collection. TFound is a modular platform and supports other annotated genomes and sets of Position Matrices, providing a simple and accessible tool to decipher multiple regulatory mechanisms through TF-cis regulatory element interaction inference, paving the way to understanding and engineering of a range of microorganisms.

## Instalation
-Install dependencies:
--Pandas (https://pandas.pydata.org/pandas-docs/stable/install.html);<br>
--Numpy (https://docs.scipy.org/doc/numpy-1.15.1/user/install.html);<br>
--MatPlotLib (https://matplotlib.org/users/installing.html);<br>
--BeautifulSoup4 (https://www.crummy.com/software/BeautifulSoup/bs4/doc/);<br>
--Numba (http://numba.pydata.org/numba-doc/0.13/install.html).<br>
-Run DbCreate.py (download the PWMs from YeTFaSCo and TF information from yeastgenome.org to create the database.)<br>
-Run InterfaceView.py
