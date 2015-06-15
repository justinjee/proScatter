# pScatter

##For visualizing pLink data from one or more experiments

by Justin Jee, with design by Katelyn McGary Shipper

![pScatter plot of interacting proteins](test.png =250x)

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

--scale

Scales both plot and output so that only the amino acids of interest are considered.

--zoom=Prot1-Prot2

Zooms in on only one subplot (Prot1 vs Prot2). As an added feature, clicking on any point in the scatter plot will print the coordinates of that point in the console.

--xkcd

Because why not.
