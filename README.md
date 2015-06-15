# pScatter

##For visualizing pLink data from one or more experiments

by Justin Jee, with design by Katelyn McGary Shipper

*Dependencies*

pScatter uses [python](https://www.python.org/downloads/), [numpy](http://www.numpy.org/), and the visualization library [matplotlib](http://matplotlib.org/)

*Basic Use*

The inputs to pScatter include:
1.   A fasta file including the amino acid sequences of all the proteins under consideration
2.   A list of amino acids of interest (ex: K or CM)
3.   One or more directories containing pLink files (in .txt format) with crosslink information (of the form ProteinA(position1)-ProteinB(position2))

The outputs include:
1.   A summary file (.txt) containing the list of all links and their frequencies
2.   A scatter plot

*Example*
python pScatter.py test.fasta K EXP1 EXP2 testout

*Features*
--fullscale

--zoom=Prot1-Prot2

--xkcd
