import pandas as pd

def get_df_orf(ref):
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

    # percent = 0.21
    # orf1a_start = None
    # orf1a_end = None
    # df_orf1a = df_orf[df_orf["stop"] > len(ref) * percent]
    # if len(df_orf1a) > 0:
    #     orf1a_start = df_orf1a.iloc[0]["start"]
    #     orf1a_end = df_orf1a.iloc[0]["stop"]

    # df_orf1ab = df_orf[(df_orf["start"] > orf1a_end - 100) & (
    #     (df_orf["stop"] - df_orf["start"]) > len(ref) * percent)]
    # if len(df_orf1ab) > 0:
    #     return orf1a_start, orf1a_end, df_orf1ab.iloc[0]["start"], df_orf1ab.iloc[0]["stop"]

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
    """Predict putative ORFs"""
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

    possible_trs = {}
    for start in start_pos:
        for stop in stop_pos:
            if stop > start and (stop - start) % 3 == 0:
                possible_trs[start] = stop
                break
        else:
            possible_trs[start] = length - 1

    return next_start, possible_trs

if __name__ == "__main__":
    seq = "ATGATCTAACAATGAATTAAGTAATCAATG"
    next_start, possible_trs = predict_ORFs(seq)
    print(''.join(f"{x:<2} " for x in range(len(seq))))
    print('  '.join(seq))
    print(''.join(f"{x:<2} " for x in next_start))
    print(possible_trs)