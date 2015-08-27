# proScatter

##For visualizing pLink data from one or more experiments

by Katelyn McGary Shipper, Justin Jee, and Ilya Shamovsky

*Dependencies*

proScatter uses [python](https://www.python.org/downloads/), [numpy](http://www.numpy.org/), and [bokeh](http://bokeh.pydata.org/en/latest/index.html), [pandas](http://pandas.pydata.org/)

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

```
     usage: proScatter.py [-h] [-a AMINOACIDS] [-z ZOOM] [-s] [-e EVALUE] [-u]
                          [-o OUTPUT] [-v]
                          fasta_file plink
     
     positional arguments:
       fasta_file            fasta file with protein sequences
       plink                 pLink output .html file
     
     optional arguments:
       -h, --help            show this help message and exit
       -a AMINOACIDS, --aminoacids AMINOACIDS
                             cross-linkable aminoacids. Defaults to Lysine (K).
       -z ZOOM, --zoom ZOOM  Prot1-Prot2 only display subplot for proteins Prot1 vs
                             Prot2
       -s, --scale           scale both plot and outputs so that only amino acids
                             of interest are considered
       -e EVALUE, --evalue EVALUE
                             e-value cutoff
       -u, --unjoin          unjoin plot axes
       -o OUTPUT, --output OUTPUT
                             output file (HTML) name
       -v, --verbose         increase output verbosity
```
