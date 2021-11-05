from functools import wraps
from time import time

def color_code(base):
    """Unix console color code for nucleotides

    Args:
        base (str): nucleotides (ACGT)

    Returns:
        str: colored nucleotides
    """
    if base == "G":
        return "\x1b[48;5;220mG\x1b[0m"
    elif base == "C":
        return "\x1b[48;5;204mC\x1b[0m"
    elif base == "A":
        return "\x1b[48;5;118mA\x1b[0m"
    elif base == "T":
        return "\x1b[48;5;81mT\x1b[0m"
    else:
        return base


def render_color_seq(seq: str):
    """Colorful base pairs

    Args:
        seq (str): base pairs

    Returns:
        str: base pairs with colors
    """
    return ''.join(map(color_code, seq))


def make_score_func(match=1, mismatch=-1, indel=-1):
    """make a scoring function

    Args:
        match (int, optional): matching score. Defaults to 1.
        mismatch (int, optional): mismatch score. Defaults to -1.
        indel (int, optional): indel score. Defaults to -1.

    Returns:
        function: scoring function
    """
    characters = ['A', 'C', 'G', 'T', 'N', 'R', 'W', 'S',
                  'Y', 'M', 'K', 'V', 'H', 'D', 'B']
    represents = {
        'A': {'A'},
        'C': {'C'},
        'G': {'G'},
        'T': {'T'},
        'N': {'A', 'C', 'G', 'T'},
        'R': {'A', 'G'},
        'W': {'A', 'T'},
        'S': {'C', 'G'},
        'Y': {'C', 'T'},
        'M': {'A', 'C'},
        'K': {'G', 'T'},
        'V': {'A', 'C', 'G'},
        'H': {'A', 'C', 'T'},
        'D': {'A', 'G', 'T'},
        'B': {'C', 'G', 'T'},
    }
    score = {}
    for c in characters:
        score[c] = {k: match if k in represents[c]
                    else mismatch for k in characters}
        score[c]['-'] = indel
    score['-'] = {k: indel for k in characters}
    score['-']['-'] = 0
    return lambda a, b: score[a][b]

def get_description(fasta):
    """get description from FASTA file

    Args:
        fasta (str): path to FASTA file

    Returns:
        str: description
    """
    with open(fasta) as ifile:
        for line in ifile:
            if line.startswith(">"):
                name_split = line.find(' ')
                return line[name_split+1:].strip()

def get_name(fasta):
    """get genome name from FASTA file

    Args:
        fasta (str): path to FASTA file

    Returns:
        str: genome name
    """
    with open(fasta) as ifile:
        for line in ifile:
            if line.startswith(">"):
                name_split = line.find(' ')
                return line[1:name_split].strip()


def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print(f'func:{f.__name__} took: {te-ts:2.4f} sec')
        return result
    return wrap