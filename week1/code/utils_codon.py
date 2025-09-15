import os
from typing import List, Tuple

def read_fasta(path: str, name: str) -> List[str]:
    data: List[str] = []
    with open(os.path.join(path, name), 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith('>'):
                continue
            data.append(line)
    print(name, len(data), len(data[0]) if data else 0)
    return data

def read_data(path: str) -> Tuple[List[str], List[str], List[str]]:
    short1: List[str] = read_fasta(path, "short_1.fasta")
    short2: List[str] = read_fasta(path, "short_2.fasta")
    long1: List[str] = read_fasta(path, "long.fasta")
    return short1, short2, long1