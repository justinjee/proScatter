#fasta
import re
import sys

def minGreaterThan(ar,n):
    m = max(ar)
    for i in ar:
        if i<m and i>n:
            m = i
    return m

#first argument: fasta filename, 2nd argument: xwalk filename, 3rd argument: a.a. of interest

f = open(sys.argv[1],'r')
#function to open and read a fasta file. remember to press i
# now read file lines one by one
#make dictionary
fd = {} #a.a. locations
ld = {} #total lengths of all proteins
name = ''
seq = ''
for line in f:
    line = line.rstrip()
    if '>' in line:
        klocs = [m.start()+1 for m in re.finditer(sys.argv[3],seq)]
        fd[name] = klocs
        ld[name] = len(seq)
        seq = ''
        name = line[1:]
    else:
        line = "".join(line.split())
        seq += line
klocs = [m.start()+1 for m in re.finditer(sys.argv[3],seq)]
fd[name] = klocs
ld[name] = len(seq)


#xwalk
f = open(sys.argv[2],'r')

xd = {}

for line in f:
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
#x2prot = {} #xlink label (A or B): 
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
                if score>maxscore:
                    maxscore = score
                    maxname = k2
                    maxoffset = i
        print k + " part "+str(part)+" is " + maxname + " at "+ str(maxoffset)
        totalscore += maxscore
        startoffset = maxoffset+ld[maxname]
        stopoffset = minGreaterThan(s1,startoffset)
        part+=1

#write new plink file

