def read_fasta(path: str, name: str):
    data = []
    filename = f"{path}/{name}"
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] != '>':
                data.append(line)
    # optional sanity print (comment out in CI if noisy)
    # print(name, len(data), len(data[0]) if data else 0)
    return data

def read_data(path: str):
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1  = read_fasta(path, "long.fasta")
    return short1, short2, long1
