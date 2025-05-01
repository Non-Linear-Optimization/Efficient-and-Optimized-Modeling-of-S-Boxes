import numpy as np
import json

def anf(s, input_bit_size):
    s = list(s)
    for i in range(input_bit_size):
        for j in range(0, 2**input_bit_size, 2*2**i):
            for o in range(2**i):
                s[j+o+2**i] ^= s[j+o]
    return s

def expand(lst, input_bit_size):
    res = set()
    for v in range(2**input_bit_size):
        if any(v & u == u for u in lst):
            res.add(v)
    return sorted(res)

def get_sbox(raw_sbox):
    # Name of the S box
    name = raw_sbox["name"]

    # Input and output size of the S box
    input_bit_size = int(raw_sbox["input"])
    output_bit_size = int(raw_sbox["output"]) 

    # S box in the form of an integer array
    s_box = [0]*(2**input_bit_size)
    for i in raw_sbox["s-box"]:
        s_box[int(i)] = int(raw_sbox["s-box"][i])
    return name, input_bit_size, output_bit_size, s_box


def gen_DPT(input_bit_size, output_bit_size, s_box):
    dppt = [set() for u in range(2**input_bit_size)]
    for v in range(2**input_bit_size):
        a = anf([int(y & v == v) for y in s_box], input_bit_size)
        for u, take in enumerate(a):
            if take:
                dppt[u].add(v)

    for u1 in range(2**input_bit_size):
        for u2 in range(u1):
            if u2 & u1 == u2: 
                dppt[u2] = dppt[u2] | dppt[u1]

    dppt = [expand(lst, input_bit_size) for lst in dppt]
    DPT = np.zeros((2**input_bit_size, 2**output_bit_size))
    DPT = DPT.astype(int)
    for u, lst in enumerate(dppt):
        for v in lst:
            DPT[u][v] = 1

    return DPT


# All S boxes are loaded into 'data'
f = open('present_sbox.json')
data = json.load(f)
f.close()

for val in data:
    # Take S box from file input
    name, input_bit_size, output_bit_size, s_box = get_sbox(val)

    # Generate DPT
    DPT = gen_DPT(input_bit_size, output_bit_size, s_box)

    print("\nDPT:")
    print(DPT)

    # possible_transitions, impossible_transitions = get_transitions(input_bit_size, output_bit_size, DPT)