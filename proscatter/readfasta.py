import re

f = open('test.fasta','r')

d = {}
name = ''
seq = ''
for line in f:
    line = line.rstrip()
    if '>' in line:
        klocs = [(m.start()+1) for m in re.finditer('K',seq)]
        d[name] = klocs
        seq = ''
        name = line[1:]
    else:
        line = "".join(line.split())
        seq += line

klocs = []
for i in range(len(seq)):
    if seq[i]=='K':
        klocs.append(i+1)
d[name] = klocs

print d
