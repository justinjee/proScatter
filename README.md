# proScatter

##For visualizing pLink data from one or more experiments

by Katelyn McGary Shipper and Justin Jee

*Dependencies*

proScatter uses [python](https://www.python.org/downloads/), [numpy](http://www.numpy.org/), and [bokeh](http://bokeh.pydata.org/en/latest/index.html)

*Basic Use*

To download, either download all 3 .py files into the same directory, or download the zipped folder (see button on the bottom right) 

The inputs to proScatter include:

1.   A fasta file including the amino acid sequences of all the proteins under consideration.
2.   A list of amino acids of interest (ex: K or CM)
3.   Up to three directories containing pLink files (in .txt format) with crosslink information (of the form ProteinA(position1)-ProteinB(position2))
     - Specifically, the file should contain these columns: 1,1	File2335 Spectrum451 scans: 7430.dta	3.30E-02	2137.126654	0.00148	0.692192	null	1.pFind	1.pFind	1	RpoC(953)-RpoC(992)
     - However, only the e value (column 6) and the crosslink (last column) are used by proScatter

Note: If converting from excel to .txt files, it is important to save it in WINDOWS TXT FORMAT. Even if you are on Mac. 
The exact features the program looks for in a given line are in the conditional "if ',' in line and not ('REVERSE' in line):"
You can modify this in loadfiles.py if your file looks different.

The outputs include:

1.   A summary file (.txt) containing the list of all links and their frequencies
2.   A scatter plot

*Example*

python proScatter.py test.fasta K EXP1 EXP2 testout

##Features

proScatter enables multiple features, for example:

python proScatter.py test.fasta K EXP1 EXP2 testout --scale --zoom=Prot1-Prot2 --evalue=0.001

Details are given below:

*--scale*

Scales both plot and output so that only the amino acids of interest are considered. Axes are in units of amino acids of interest (ex: 1st lysine, 2nd lysine, etc)

*--zoom=Prot1-Prot2*

Zooms in on only one subplot (Prot1 vs Prot2). As an added feature, clicking on any point in the scatter plot will print the coordinates of that point in the console.

*--evalue=#*

Considers only links with a score below a certain number #. Scores are expected to be in the 5th column
