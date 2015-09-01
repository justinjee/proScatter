# proScatter

##For visualizing pLink data from one or more experiments

by Katelyn McGary Shipper, Justin Jee, and Ilya Shamovsky

*Dependencies*

proScatter uses [python](https://www.python.org/downloads/), [numpy](http://www.numpy.org/), [bokeh](http://bokeh.pydata.org/en/latest/index.html), and [pandas](http://pandas.pydata.org/)

*Basic Use*

To download, either download all 3 .py files into the same directory, or download the zipped folder (see button on the bottom right) 

The inputs to proScatter include:

1.   A fasta file including the amino acid sequences of all the proteins under consideration.
2.   A list of amino acids of interest (ex: K or CM)
3.   pLink output in `.html` format

The outputs include:

1.   A summary file (.txt) containing the list of all links and their frequencies
2.   A scatter plot

*Example*

`./proScatter.py examples/test.fasta examples/IdProTable_combine.html`

## Installation

Clone the repo, `cd` into the cloned dir and run:

`python setup.py install`


##Features

proScatter enables multiple features, for example:

`python proScatter.py test.fasta K test.html --scale --zoom=Prot1-Prot2 --evalue=0.001`

Details are given below:

*--scale*

Scales both plot and output so that only the amino acids of interest are considered. Axes are in units of amino acids of interest (ex: 1st lysine, 2nd lysine, etc)

*--zoom=Prot1-Prot2*

Zooms in on only one subplot (Prot1 vs Prot2). As an added feature, clicking on any point in the scatter plot will print the coordinates of that point in the console.

*--evalue=#*

Considers only links with a score below a certain number #. Scores are expected to be in the 5th column

