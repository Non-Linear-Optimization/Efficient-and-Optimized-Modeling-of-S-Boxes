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
            alpha = bin(i)[2:].zfill(input_bit_size)
            beta = bin(j)[2:].zfill(input_bit_size)
            point = [int(bit) for bit in alpha + beta]
            if DPT[i][j] != 0:
                possible_transitions.append(point)
            else:
                impossible_transitions.append(point)
    
    return possible_transitions, impossible_transitions


def gen_inequalities(possible_transitions):
    # Vertex representation of the possible transitions of the DPT
    P = Polyhedron(vertices=possible_transitions)

    # Hyperplane representation i.e., inequalities
    inequalities = list(P.Hrepresentation())

    return inequalities


def reduce_inequalities_first(impossible_transitions, inequalities):
    final_ineqs = list()
    # Greedy Approach
    while(1):
        max_count = 0
        # For each inequality count the number of impossible transitions it can remove
        inequality_lists = {q: list() for q in inequalities}
        for q in inequalities:
            count = 0
            for point in impossible_transitions:
                sum = q[0]
                for i in range(len(point)):
                    sum += (q[i+1]*point[i])
                if sum<0:
                    inequality_lists[q].append(point)
                    count += 1
            if count > max_count:
                max_count = count
        if(max_count==0):
            break

        # Obtaining the inequalities which can remove the maximum number of impossible transitions
        max_ineqaulities_list = []
        for q in inequalities:
            if len(inequality_lists[q])==max_count:
                max_ineqaulities_list.append(q)
                # max_ineqaulity = q
        
        # Select one among these max inequalities
        # First max
        q = max_ineqaulities_list[0]

        # # Mid max
        # q = max_ineqaulities_list[(len(max_ineqaulities_list)-1)//2]

        # # Last max
        # q = max_ineqaulities_list[(len(max_ineqaulities_list)-1)]

        # Random max
        # x = py_rand.randint(0,len(max_ineqaulities_list)-1)
        # q = max_ineqaulities_list[x]
        max_ineqaulity = q

        # Removing the max inequality from the list of inequalities and adding this to our final result
        inequalities.remove(max_ineqaulity)
        final_ineqs.append(max_ineqaulity)

        # Impossible transitions removed by this inequality are to be removed from the set of impossible transitions
        for p in inequality_lists[max_ineqaulity]:
            impossible_transitions.remove(p)
    return final_ineqs


def reduce_inequalities_last(impossible_transitions, inequalities):
    final_ineqs = list()
    # Greedy Approach
    while(1):
        max_count = 0
        # For each inequality count the number of impossible transitions it can remove
        inequality_lists = {q: list() for q in inequalities}
        for q in inequalities:
            count = 0
            for point in impossible_transitions:
                sum = q[0]
                for i in range(len(point)):
                    sum += (q[i+1]*point[i])
                if sum<0:
                    inequality_lists[q].append(point)
                    count += 1
            if count > max_count:
                max_count = count
        if(max_count==0):
            break

        # Obtaining the inequalities which can remove the maximum number of impossible transitions
        max_ineqaulities_list = []
        for q in inequalities:
            if len(inequality_lists[q])==max_count:
                max_ineqaulities_list.append(q)
                # max_ineqaulity = q
        
        # Select one among these max inequalities
        # # First max
        # q = max_ineqaulities_list[0]

        # # Mid max
        # q = max_ineqaulities_list[(len(max_ineqaulities_list)-1)//2]

        # Last max
        q = max_ineqaulities_list[(len(max_ineqaulities_list)-1)]

        # # Random max
        # x = py_rand.randint(0,len(max_ineqaulities_list)-1)
        # q = max_ineqaulities_list[x]
        max_ineqaulity = q

        # Removing the max inequality from the list of inequalities and adding this to our final result
        inequalities.remove(max_ineqaulity)
        final_ineqs.append(max_ineqaulity)

        # Impossible transitions removed by this inequality are to be removed from the set of impossible transitions
        for p in inequality_lists[max_ineqaulity]:
            impossible_transitions.remove(p)
    return final_ineqs


def reduce_inequalities_mid(impossible_transitions, inequalities):
    final_ineqs = list()
    # Greedy Approach
    while(1):
        max_count = 0
        # For each inequality count the number of impossible transitions it can remove
        inequality_lists = {q: list() for q in inequalities}
        for q in inequalities:
            count = 0
            for point in impossible_transitions:
                sum = q[0]
                for i in range(len(point)):
                    sum += (q[i+1]*point[i])
                if sum<0:
                    inequality_lists[q].append(point)
                    count += 1
            if count > max_count:
                max_count = count
        if(max_count==0):
            break

        # Obtaining the inequalities which can remove the maximum number of impossible transitions
        max_ineqaulities_list = []
        for q in inequalities:
            if len(inequality_lists[q])==max_count:
                max_ineqaulities_list.append(q)
                # max_ineqaulity = q
        
        # Select one among these max inequalities
        # # First max
        # q = max_ineqaulities_list[0]

        # Mid max
        q = max_ineqaulities_list[(len(max_ineqaulities_list)-1)//2]

        # # Last max
        # q = max_ineqaulities_list[(len(max_ineqaulities_list)-1)//2]

        # Random max
        # x = py_rand.randint(0,len(max_ineqaulities_list)-1)
        # q = max_ineqaulities_list[x]
        max_ineqaulity = q

        # Removing the max inequality from the list of inequalities and adding this to our final result
        inequalities.remove(max_ineqaulity)
        final_ineqs.append(max_ineqaulity)

        # Impossible transitions removed by this inequality are to be removed from the set of impossible transitions
        for p in inequality_lists[max_ineqaulity]:
            impossible_transitions.remove(p)
    return final_ineqs


def reduce_inequalities_rand(impossible_transitions, inequalities):
    final_ineqs = list()
    rand_list = ""
    # Greedy Approach
    while(1):
        max_count = 0
        # For each inequality count the number of impossible transitions it can remove
        inequality_lists = {q: list() for q in inequalities}
        for q in inequalities:
            count = 0
            for point in impossible_transitions:
                sum = q[0]
                for i in range(len(point)):
                    sum += (q[i+1]*point[i])
                if sum<0:
                    inequality_lists[q].append(point)
                    count += 1
            if count > max_count:
                max_count = count
        if(max_count==0):
            break

        # Obtaining the inequalities which can remove the maximum number of impossible transitions
        max_ineqaulities_list = []
        for q in inequalities:
            if len(inequality_lists[q])==max_count:
                max_ineqaulities_list.append(q)
                # max_ineqaulity = q
        
        # Select one among these max inequalities
        # # First max
        # q = max_ineqaulities_list[0]

        # # Mid max
        # q = max_ineqaulities_list[(len(max_ineqaulities_list)-1)//2]

        # # Last max
        # q = max_ineqaulities_list[(len(max_ineqaulities_list)-1)//2]

        # Random max
        x = py_rand.randint(0,len(max_ineqaulities_list)-1)
        rand_list += f" {x}:{len(max_ineqaulities_list)} "
        q = max_ineqaulities_list[x]
        max_ineqaulity = q

        # Removing the max inequality from the list of inequalities and adding this to our final result
        inequalities.remove(max_ineqaulity)
        final_ineqs.append(max_ineqaulity)

        # Impossible transitions removed by this inequality are to be removed from the set of impossible transitions
        for p in inequality_lists[max_ineqaulity]:
            impossible_transitions.remove(p)
    return final_ineqs, rand_list


def check_inequalities(possible_transitions, impossible_transitions, final_ineqs):
    for q in final_ineqs:
        for point in possible_transitions:
            sum = q[0]
            for i in range(len(point)):
                sum += (q[i+1]*point[i])
            if sum<0:
                return False
        count = 0
        for point in impossible_transitions:
            sum = q[0]
            for i in range(len(point)):
                sum += (q[i+1]*point[i])
            if sum>=0:
                count += 1
        if (count == len(impossible_transitions)):
            return False
    return True


if __name__ == "__main__":

    # All S boxes are loaded into 'data'
    f = open('../SBOXES/5_bit_sboxes.json')
    data = json.load(f)
    f.close()

    
        # for q in final_ineqs_first:
        #     f1.write(str(list(q))+"\n")
        # f1.close()

    for val in data:
        # Take S box from file input
        name, input_bit_size, output_bit_size, s_box = get_sbox(val)

        # Generate DPT
        DPT = gen_DPT(input_bit_size, output_bit_size, s_box)

        # print(DPT)

        # Get possible and impossible transitions
        possible_transitions, impossible_transitions = get_transitions(input_bit_size, output_bit_size, DPT)

        impossible_transitions_1 = impossible_transitions.copy()
        # impossible_transitions_2 = impossible_transitions.copy()
        # impossible_transitions_3 = impossible_transitions.copy()
        # impossible_transitions_5 = impossible_transitions.copy()

        # Generate Inequalities using sage
        inequalities = gen_inequalities(possible_transitions)

        # Reduce the inequalities
        final_ineqs_first = reduce_inequalities_first(impossible_transitions_1, inequalities)
        # final_ineqs_last = reduce_inequalities_last(impossible_transitions_2, inequalities)
        # final_ineqs_mid = reduce_inequalities_mid(impossible_transitions_3, inequalities)
        # final_ineqs_rand, rand_list = reduce_inequalities_rand(impossible_transitions_5, inequalities)

        f2 = open(f"DPT_Greedy_5-bit_sboxes.txt","a")
        f2.write(name + " - " + str(len(final_ineqs_first)) + "\n")
        f2.close()

        f1 = open(f"5-bit_sboxes/{name}_Greedy_DPT.txt","a")
        for q in final_ineqs_first:
            f1.write(str(list(q))+"\n")
        f1.close()

        
        # f1 = open(f"{name}_last.txt","a")
        # for q in final_ineqs_last:
        #     f1.write(str(list(q))+"\n")
        # f1.close()
        # f1 = open(f"{name}_mid.txt","a")
        # for q in final_ineqs_mid:
        #     f1.write(str(list(q))+"\n")
        # f1.close()
        # f1 = open(f"Random/{name}_rand.txt","a")
        # for q in final_ineqs_last:
        #     f1.write(str(list(q))+"\n")
        # f1.write(f"Order chosen for this result- index of chosen inequality:number of max inequalities\n")
        # f1.write(f"{str(rand_list)}")
        # f1.close()
    