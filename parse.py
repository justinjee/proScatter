import re, HTMLParser

RE = re.compile(r'(?P<segid1>\w+)[(](?P<resid1>\d+)[)]-(?P<segid2>\w+)[(](?P<resid2>\d+)[)]')        

class PLinkParser(HTMLParser.HTMLParser):

    def __init__(self):
        #super(TxtParser, self).__init__()
        HTMLParser.HTMLParser.__init__(self)
        self.buf = []
        self._line = ''
        self._row = False
        self._td = False
        self._td_count = 0
        self._grab = False
        
    def handle_starttag(self, tag, attrs):
        ntag = tag.upper()
        if ntag == 'TR':
            self._row = True
        elif ntag == 'TD':
            self._td = True
            self._td_count += 1
        
    def handle_endtag(self, tag):
        ntag = tag.upper()
        if ntag == 'TR':
            self._row = False
            self._td_count = 0
            if self._line:
                self.buf.append(self._line)
            self._line = ''
        elif ntag == 'TD':
            self._td = False
        
    def handle_data(self, data):
        if self._td and self._td_count == 1 and data.strip().isdigit():
            self._grab = True
        elif self._td_count == 2 and self._grab:
            self._line += data.strip()
            self._grab = False
        
def extract_xl_coord(line):
    m = RE.search(line)
    return {
        'prot1': m.group('segid1'),
        'resid1': m.group('resid1'),
        'prot2': m.group('segid2'),
        'resid2': m.group('resid2'),
        }
    
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
