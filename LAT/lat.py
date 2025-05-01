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


def gen_LAT(input_bit_size, output_bit_size, s_box):
    # Generate LAT from the obtained S-box
    LAT = np.zeros((2**input_bit_size, 2**output_bit_size))
    LAT = LAT.astype(int)
    for alpha in range(2**input_bit_size):
        for beta in range(2**output_bit_size):
            count = 0
            for x in range(2**input_bit_size):
                input_parity = bin(alpha & x).count('1') % 2
                output_parity = bin(beta & s_box[x]).count('1') % 2
                if input_parity == output_parity:
                    count += 1
            LAT[alpha][beta] = count - (2**(input_bit_size-1))
    return LAT



# All S boxes are loaded into 'data'
f = open('present_sbox.json')
data = json.load(f)
f.close()

for val in data:
    # Take S box from file input
    name, input_bit_size, output_bit_size, s_box = get_sbox(val)

    # Generate DPT
    LAT = gen_LAT(input_bit_size, output_bit_size, s_box)

    print("\nLAT:")
    print(LAT)