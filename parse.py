import re, HTMLParser
import numpy as np

RE = re.compile(r'(?P<segid1>\w+)[(](?P<resid1>\d+)[)]-(?P<segid2>\w+)[(](?P<resid2>\d+)[)]')        

HEADER_SUM = [
    ('Order', int),
    ('ProteinAC', str),
    ('MW', float),
    ('pI', float),
    ('UniquePepNum', int),
    ('Coverage(%)', float),
    ('SpecNum', int),
    ('Non-ModifiedSpecNum', int),
    ('ModifiedSpecNum', int),
    ('UniqueModifiedPepNum', int),
    ('Description', str)
]


HEADER_DETAILS = [
    #'',
    ('Order', int),
    ('Spectrum', str),
    ('Sequence', str),
    ('Score', float),
    ('Calc_M', float),
    ('Delta_M', float),
    ('ppm', float),
    ('Modification', str),
    ('Sample', str),
    ('Engine', str),
    ('MatchedIons', int),
    ('MissCleaveNum', int),
    ('Rank', int),
    ('Proteins', str)
]

def check_lists(lst1, lst2):
    if len(lst1) != len(lst2):
        return False
    return len([x for x,y in zip(lst1, lst2) if x != y]) == 0


class PLinkParser(HTMLParser.HTMLParser):

    def __init__(self):
        HTMLParser.HTMLParser.__init__(self)
        self.sum_records = []
        self.details_records = []
        self._current = None
        self._current_header = None
        self.buf = []
        self._td = False
        
    def handle_starttag(self, tag, attrs):
        ntag = tag.upper()
        if ntag == 'TD':
            self._td = True
                
    def _add_record(self):
        record = {}
        cleaned = [x for x in self.buf if x != '']
        for key,value in zip(self._current_header, cleaned):
            try:
                record[key[0]] = key[1](value)
            except ValueError:
                record[key[0]] = np.nan
        self._current.append(record)
        self.buf = []
        
    def handle_endtag(self, tag):
        ntag = tag.upper()
        if ntag == 'TR':
            if check_lists([x[0] for x in HEADER_SUM], self.buf):
                self._current, self._current_header = self.sum_records, HEADER_SUM
            elif check_lists([x[0] for x in HEADER_DETAILS], self.buf):
                self._current, self._current_header = self.details_records, HEADER_DETAILS
            else:
                self._add_record()
            self.buf = []
        elif ntag == 'TD':
            self._td = False
            
    def handle_data(self, data):
        if self._td and data:
            self.buf.append(data.strip())
    
    
    
def is_empty(line):
    """Returns True empty lines and lines consisting only of whitespace."""
    return (not line) or line.isspace()

def is_fasta_label(line):
    return line.strip().startswith('>')
    

def LabeledRecordFinder(is_label_line, ignore=is_empty):
    """Returns function that returns successive labeled records from file.
    Includes label line in return value. Returns list of relevant lines.
    Skips over any lines for which ignore(line) evaluates True (default is
    to skip empty lines).
    """
    def parser(lines):
        with open(lines, 'r') as lines:
            curr = []
            for l in lines:
                line = l.strip()
                if ignore(line):
                    continue
                # if we find the label, return the previous record
                if is_label_line(line):
                    if curr:
                        yield curr
                        curr = []
                curr.append(line)
            # don't forget to return the last record in the file
            if curr:
                yield curr
    return parser

    
FastaFinder = LabeledRecordFinder(is_fasta_label)


def parse_fasta(infile, finder=FastaFinder):
    """Generator of labels and sequences from a fasta file.
    """

    for rec in finder(infile):
        # first line must be a label line
        if not rec[0].startswith('>'):
            raise ValueError("Found Fasta record without label line: %s" % rec)
        # record must have at least one sequence
        if len(rec) < 2:
            raise ValueError("Found label line without sequences: %s" % rec)

        # remove the label character from the beginning of the label
        label = rec[0][1:].strip()
        seq = ''.join(rec[1:])

        yield label, seq
