import numpy as np

def one_hot_encode(seq):
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    one_hot_encodings = np.eye(4)
    seq_to_mappings = [mapping[base] for base in seq]
    return(one_hot_encodings[seq_to_mappings])
