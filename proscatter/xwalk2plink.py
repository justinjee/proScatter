#fasta
import re
import sys
import parse
from collections import defaultdict

#first argument: fasta filename, 2nd argument: xwalk filename, 3rd argument: a.a. of interest
#EXAMPLE: python readxwalk2.py ECRNINVIVO1.fasta uvrdrnap.txt K > output.txt

parser = argparse.ArgumentParser(description='''
    pScatter
    Interaction visualizer for pLink XLMS data, by Justin Jee (with design by Katelyn McGary Shipper)

    '''
    )

parser.add_argument('fasta_file', type=str, help='fasta file with protein sequences')
parser.add_argument('xwalk_file', type=str, help='XWalk output file')
parser.add_argument('-a', '--aminoacid', default='K', help='cross-linkable aminoacids. Defaults to Lysine (K).')
args = parser.parse_args()

fd = {}
ld = {}
with open(args.fasta_file, 'r') as fi:
    for name,seq in parse.parse_fasta(fi):
        fd[name] = [m.start() + 1 for m in re.finditer(args.aminoacid, seq)]
        ld[name] = len(seq)

#xwalk
# Line
# 1 uvrdnap_ecmodel.pdb LYS-486-A-CB    LYS-496-B-CB    652 10.7    11.5    -   -   -

with open(args.xwalk_file, 'r') as fi:

    xd = {}

    for line in fi:
        temparray = line.split()
        temp1 = temparray[2].split('-')
        temp2 = temparray[3].split('-')

# instantiate     
        if not temp1[2] in xd:
            xd[temp1[2]] = set([])
        if not temp2[2] in xd:
            xd[temp2[2]] = set([])    
        xd[temp1[2]].add(int(temp1[1]))
        xd[temp2[2]].add(int(temp2[1]))
    
#compare
    x2prot = defaultdict(list) #xlink label (A or B): sorted array of (end coordiate, protein name) 
    for k in xd:
        s1 = xd[k]
        totalscore = 0
        maxname = ''
        part = 1
        startoffset=0
        stopoffset=0
        while totalscore < len(s1):
            maxscore = 0
            maxoffset = 0
            for i in range(startoffset,stopoffset+1):
                si = set([j-i for j in s1])
                for k2 in fd:
                    score = len(si.intersection(fd[k2]))
                    if score >  maxscore:
                        maxscore = score
                        maxname = k2
                        maxoffset = i
            #print k + " part "+str(part)+" is " + maxname + " at "+ str(maxoffset)
            x2prot[k].append((maxoffset,maxname))
            totalscore += maxscore
            startoffset = maxoffset+ld[maxname]
            stopoffset = minGreaterThan(s1,startoffset)
            part += 1

#write new plink file

def translate(name, loc, x2prot):
    ar = x2prot[name]
    for i,_ in enumerate(ar):
        if (i >= length - 1) or (loc < ar[i+1][0]):
            newname = ar[i][1]
            newloc = loc - ar[i][0]
    return '{0}({1})'.format(newname, newloc)

with open(args.xwalk_file, 'r') as fi:
    for line in fi:
        temparray = line.split()
        temp1 = temparray[2].split('-')
        temp2 = temparray[3].split('-')
        (loc1, name1) = temp1[1:3]
        (loc2, name2) = temp2[1:3]
        print ', # # # # 0 # # # # # # # '+translate(name1,int(loc1),x2prot)+'-'+translate(name2,int(loc2),x2prot)
