import numpy as np
import json
import time
import gurobipy as gp
from gurobipy import GRB
from math import gcd
from functools import reduce
from itertools import combinations

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

def get_transitions(input_bit_size, output_bit_size, DPT):
    # Get possible and impossible transitions from the DPT
    possible_transitions = []
    impossible_transitions = []    
    
    for i in range(2**input_bit_size):
        for j in range(2**input_bit_size):
            x = bin(i)[2:].zfill(input_bit_size)
            y = bin(j)[2:].zfill(input_bit_size)
            point = int(x + y, 2)
            if DPT[i][j] != 0:
                possible_transitions.append(point)
            else:
                impossible_transitions.append(point)
    
    return possible_transitions, impossible_transitions


def gen_function(possible_transitions, impossible_transitions, n, sbox_name, a_bound, b_bound):
    B = set(impossible_transitions)
    p = ''
    Final_inequalities = list()
    while B:
        # Create Model
        M = gp.Model()

        # Variables
        a_vars = M.addVars(n, vtype=GRB.INTEGER, lb=-a_bound, ub=a_bound, name="a")
        b = M.addVar(vtype=GRB.INTEGER, lb=0, ub=b_bound, name="b")
        y_vars = M.addVars(B, vtype=GRB.BINARY, name="y")
        # z_vars = M.addVars(border, vtype=GRB.BINARY, name="z")

        # Constraints
        for v in possible_transitions:
            bin_v = format(v, f'0{n}b')
            f_v = sum(a_vars[i] * int(bin_v[i]) for i in range(n)) + b
            M.addConstr(f_v >= 0, name=f"constraint_1_{v}")
        for v in B:
            bin_v = format(v, f'0{n}b')
            f_v = sum(a_vars[i] * int(bin_v[i]) for i in range(n)) + b
            M.addConstr(f_v - ((n*a_bound + b_bound + 1)*(1-y_vars[v])) <= -1, name=f"constraint_2_{v}")



        # Objective
        M.setObjective(sum(y_vars[v] for v in B), GRB.MAXIMIZE)

        # Optimize
        M.optimize()

        Result = [int(a_vars[i].X) for i in range(n)] + [int(b.X)]
        Final_inequalities.append(Result)
        if M.status == GRB.OPTIMAL:
            for v in list(B):
                if int(y_vars[v].X) == 1:
                    B.discard(v)

        # Dispose
        M.dispose()

    return Final_inequalities

if __name__ == "__main__":

    f = open('../SBOXES/4_bit_sboxes.json')
    # Take S boxes from file input
    data = json.load(f)
    f.close()

    for val in data:
        # Load S box from data
        name, input_bit_size, output_bit_size, s_box = get_sbox(val)

        # Generate DPT
        DPT = gen_DPT(input_bit_size, output_bit_size, s_box)

        # Get the possible and impossible transitions
        possible_transitions, impossible_transitions = get_transitions(input_bit_size, output_bit_size, DPT) 

        # Intialize variables for this method
        n = input_bit_size+output_bit_size

        # Select the bounds
        a_bound = 2**(input_bit_size-1)
        b_bound = 2**(input_bit_size)

        # Final Inequalities
        Final_inequalities = gen_function(possible_transitions, impossible_transitions, n, name, a_bound, b_bound)
        
        file1 = open(f"4-bit_sboxes/{name}.txt", "a")
        for q in Final_inequalities:
            file1.write(f"{str(list(q))}\n")
        file1.close()

        f2 = open(f"DPT_DIG_4-bit_sboxes.txt","a")
        f2.write(name + " - " + str(len(Final_inequalities)) + "\n")
        f2.close()
