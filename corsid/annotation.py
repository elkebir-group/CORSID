import pandas as pd


def parse_attributes(attr: str):
    l = []
    for pair in attr.split(";"):
        if len(pair) > 0:
            k, v = pair.split('=')
        l.append([k.lower(), v])
    return dict(l)


def get_annotation_region(annotation_file: str):
    annot = pd.read_csv(annotation_file,
                        sep='\t',
                        header=None,
                        names=["id", "source", "feature", "start", "end",
                               "score", "strand", "phase", "attributes"],
                        comment='#')
    regions = {}
    for _, row in annot[annot["feature"] == "CDS"].iterrows():
        infos = parse_attributes(row["attributes"])
        if "gene" in infos:
            name = infos["gene"].lower()
        elif "product" in infos:
            name = infos["product"].lower()
        elif "pseudo" in infos and infos["pseudo"] == "true":
            continue
        else:
            print("No label: ", row)
            exit(1)

        known_name_of_1ab = {"1a", "1b", "1ab",
                             "replicase", "polymerase", "polyprotein", "pol"}
        if any(n in name.lower() for n in known_name_of_1ab) or name.lower() == '1' or name.lower() == 'orf1':
            continue

        if name in convert:
            name = convert[name]
        if name.startswith("hypothetical protein"):
            name = name[20:].strip() + '?'
        if name.startswith("non-structural protein"):
            name = "ns" + name[22:].strip()
        if name.startswith("nonstructural protein"):
            name = "ns" + name[21:].strip()
        if name.startswith("putative"):
            name = name[8:].strip() + '??'
        if name.endswith("protein"):
            name = name[:-7].strip()
        if name.startswith("protein"):
            name = name[7:].strip()
        if name.startswith("orf"):
            name = name[3:].strip()
        if name in regions:
            # take the union region if already overlapping a 
            a = regions[name][0]
            b = regions[name][1]
            c = row["start"] - 1
            d = row["end"] - 1
            if (
                (a <= d and c <= b and
                (min(b, d) - max(a, c)) / (max(b, d) - min(a, c)) >= 0.1)
            ):
                regions[name] = (min(a, c), max(b, d))
            else:
                while name in regions:
                    name += '#'
                regions[name] = (row["start"] - 1, row["end"] - 1)
        else:
            regions[name] = (row["start"] - 1, row["end"] - 1)
    return regions

convert = {
    '3A protein': '3a',
    '3a protein': '3a',
    '5a protein#2': '5a#2',
    '7A': '7a',
    '7B': '7b',
    'Membrane protein': 'M',
    'Spike': "S",
    'Spike (S) glycoprotein': "S",
    'e': "E",
    'envelop protein': "E",
    'envelop protein (E)': "E",
    'envelope': "E",
    'envelope E protein': "E",
    'envelope protein': "E",
    'hemagglutinin esterase': "HE",
    'hemagglutinin-esterase': "HE",
    'hemagglutinin-esterase glycoprotein': "HE",
    'hemagglutinin-esterase glycoprotein (HE)': "HE",
    'hemagglutinin-esterase protein': "HE",
    'hypothetical protein#2': "?#2",
    'm': "M",
    'membrance glycoprotein': "M",
    'membrane': "M",
    'membrane M protein': "M",
    'membrane glycoprotein': "M",
    'membrane glycoprotein M': "M",
    'membrane protein': "M",
    'membrane protein (M)': "M",
    'n': "N",
    'nonstructural protein#2': "NS#2",
    'nucleocapsid': "N",
    'nucleocapsid N protein': "N",
    'nucleocapsid phosphoprotein': "N",
    'nucleocapsid protein': "N",
    'nucleocapsid protein (N)': "N",
    'nucleopcapsid N protein': "N",
    'nucleoprotein': "N",
    'nucloecapsid protein': "N",
    'putative internal (I) ORF': "I?",
    's': "S",
    'small envelope E protein': "small E",
    'small envelope protein': "small E",
    'small membrane protein': "small M",
    'spike': "S",
    'spike S glycoprotein': "S",
    'spike glycoprotein': "S",
    'spike glycoprotein (S)': "S",
    'spike protein': "S",
    'spike surface glycoprotein': "S"
}
