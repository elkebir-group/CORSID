import pandas as pd

def get_df_orf(ref):
    """Get as dataframe of all possible ORFs

    Args:
        ref (str): reference genome

    Returns:
        pd.DataFrame: dataframe of all possible ORFs
    """
    possible_orf = []
    start_pos = []
    start_pos.append(ref.find('ATG'))
    while start_pos[-1] != -1:
        start_pos.append(ref.find('ATG', start_pos[-1]+1))

    stop_codons = ["TAA", "TAG", "TGA"]
    stop_pos = []
    for i in range(len(ref)-3):
        if ref[i:i+3] in stop_codons:
            stop_pos.append(i)

    possible_orf = []
    for start in start_pos:
        for stop in stop_pos:
            if stop > start and (stop - start) % 3 == 0:
                possible_orf.append((start, stop))
                break

    df_orf = pd.DataFrame(possible_orf, columns=["start", "stop"])
    df_orf["length"] = df_orf["stop"] - df_orf["start"]
    df_orf.sort_values("start", inplace=True)
    return df_orf.iloc[1:, :]

def guess_orf1ab(ref):
    """Heuristic for finding ORF1ab

    Args:
        ref (str): reference genome

    Returns:
        tuple: start of 1a, end of 1a, end of 1b
    """
    df_orf = get_df_orf(ref)

    length = len(ref)
    founds = []
    stop_codons = ["TAA", "TAG", "TGA"]
    for idx, row in df_orf.iterrows():
        found = []
        for i in range(row["stop"]-1, length-3, 3):
            if ref[i:i+3] in stop_codons:
                found = (row["start"], row["stop"], i)
                break
        else:
            break
        for i in range(row["stop"]-2, length-3, 3):
            if ref[i:i+3] in stop_codons:
                if i - row["start"] > found[2] - found[0]:
                    found = (row["start"], row["stop"], i)
                break
        founds.append(found)
    max_idx = 0
    for i, f in enumerate(founds):
        if f[2] - f[0] > founds[max_idx][2] - founds[max_idx][0]:
            max_idx = i
    return founds[max_idx][0], founds[max_idx][1], founds[max_idx][2]


def predict_ORFs(seq: str):
    """Predict putative ORFs

    Args:
        seq (str): genome

    Returns:
        tuple: the list of next start codons, and corresponding stop codon for each start codon
    """
    length = len(seq)
    start_pos = []
    start_pos.append(seq.find('ATG'))
    while True:
        next_start = seq.find('ATG', start_pos[-1]+1)
        if next_start != -1:
            start_pos.append(next_start)
        else:
            break
    next_start = [None] * length
    pointer = 0
    for i in range(length):
        if i <= start_pos[pointer]+2:
            next_start[i] = start_pos[pointer]
        else:
            pointer += 1
            if pointer < len(start_pos):
                next_start[i] = start_pos[pointer]
            else:
                break

    stop_codons = ["TAA", "TAG", "TGA"]
    stop_pos = []
    for i in range(length-3):
        if seq[i:i+3] in stop_codons:
            stop_pos.append(i)

    next_stop = {}
    for start in start_pos:
        for stop in stop_pos:
            if stop > start and (stop - start) % 3 == 0:
                next_stop[start] = stop
                break
        else:
            next_stop[start] = length - 1

    return next_start, next_stop
