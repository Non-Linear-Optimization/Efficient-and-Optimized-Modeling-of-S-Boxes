import numpy as np
import json
import time
import random as py_rand
from sage.all import *

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


def gen_BCT(input_bit_size, output_bit_size, s_box):
    # Generate BCT from the obtained S-box
    BCT = np.zeros((2**input_bit_size, 2**output_bit_size), dtype=int)
    
    for alpha in range(2**input_bit_size):  # Input difference α
        for beta in range(2**output_bit_size):  # Output difference β
            count = 0
            for x in range(2**input_bit_size):
                term_1 = s_box.index(s_box[x] ^ beta)
                term_2 = s_box.index(s_box[x ^ alpha] ^ beta)
                
                if (term_1 ^ term_2) == alpha:  # Boomerang condition
                    count += 1
            BCT[alpha][beta] = count
    
    return BCT

def get_transitions(input_bit_size, output_bit_size, BCT):
    # Get possible and impossible transitions from the BCT
    possible_transitions = []
    impossible_transitions = []    
    
    for i in range(2**input_bit_size):
        for j in range(2**input_bit_size):
            alpha = bin(i)[2:].zfill(input_bit_size)
            beta = bin(j)[2:].zfill(input_bit_size)
            point = [int(bit) for bit in alpha + beta]
            if BCT[i][j] != 0:
                possible_transitions.append(point)
            else:
                impossible_transitions.append(point)
    
    return possible_transitions, impossible_transitions


if __name__ == "__main__":

    # All S boxes are loaded into 'data'
    f = open('present_sbox.json')
    data = json.load(f)
    f.close()

    for val in data:
        # Take S box from file input
        name, input_bit_size, output_bit_size, s_box = get_sbox(val)

        # Generate BCT
        BCT = gen_BCT(input_bit_size, output_bit_size, s_box)

        print(BCT)

        # Get possible and impossible transitions
        possible_transitions, impossible_transitions = get_bct_transitions(input_bit_size, output_bit_size, BCT)
