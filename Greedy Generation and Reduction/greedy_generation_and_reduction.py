import numpy as np
import json
import time
import gurobipy as gp
from gurobipy import GRB

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

def gen_DDT(input_bit_size, output_bit_size, s_box):
    # Generate DDT from the obtained S box
    DDT = np.zeros((2**input_bit_size,2**output_bit_size), dtype=int)
    for p1 in range(2**input_bit_size):
        for p2 in range(2**output_bit_size):
            XOR_IN = p1 ^ p2
            XOR_OUT = s_box[p1] ^ s_box[p2]
            DDT[XOR_IN][XOR_OUT] += 1
    return DDT

def get_transitions(input_bit_size, output_bit_size, DDT):
    # Get possible and impossible transitions from the DDT
    possible_transitions = set()
    impossible_transitions = set()    
    for i in range(2**input_bit_size):
        for j in range(2**output_bit_size):
            x = bin(i)[2:].zfill(input_bit_size)
            y = bin(j)[2:].zfill(output_bit_size)
            point = int(x + y, 2)
            if DDT[i][j] > 0:
                possible_transitions.add(point)
            else:
                impossible_transitions.add(point)
    return possible_transitions,impossible_transitions


def gen_functions(possible_transitions, impossible_transitions, n, sbox_name):
    all_ineqs = list()
    B = set(impossible_transitions)
    removed_set = set()

    # Create Model
    M = gp.Model()

    a_bound = n
    b_bound = 2**(n/2)

    # Variables
    a_vars = M.addVars(n, vtype=GRB.INTEGER, lb=-a_bound , ub=a_bound, name="a")
    b = M.addVar(vtype=GRB.INTEGER, lb=0, ub=b_bound , name="b")
    y_vars = M.addVars(B, vtype=GRB.BINARY, name="y")

    # Constraints
    for v in possible_transitions:
        bin_v = format(v, f'0{n}b')
        f_v = sum(a_vars[i] * int(bin_v[i]) for i in range(n)) + b
        M.addConstr(f_v >= 0, name=f"constraint_1_{v}")

    for v in B:
        bin_v = format(v, f'0{n}b')
        f_v = sum(a_vars[i] * int(bin_v[i]) for i in range(n)) + b
        M.addConstr(f_v - (((n*a_bound)+ b_bound + 1)*(1-y_vars[v])) <= -1, name=f"constraint_2_{v}")

    # Objective
    M.setObjective(sum(y_vars[v] for v in B), GRB.MAXIMIZE)
    count = 0
    
    while True :
    # while len(removed_set) != len(impossible_transitions):
        N = set()
        s = ''
        
        # Optimize
        M.optimize()
        if M.status == GRB.INFEASIBLE :
            break

        q = [int(a_vars[i].X) for i in range(n)] + [int(b.X)]
        all_ineqs.append(q)
        count += 1

        # to_remove = []
        if M.status == GRB.OPTIMAL:
            for v in list(B):
                if int(y_vars[v].X) == 1:
                    s = s + str(v) + " "
                    N.add(v)
                    removed_set.add(v)

        M.addConstr(sum(y_vars[v] for v in N) <= (len(N)-1), name=f"constraint_3_{count}")
        M.addConstr(sum(y_vars[v] for v in (B - N)) >= 1, name=f"constraint_4_{count}")

        f = open(f"{sbox_name}_All_Candid_ineqs.txt","a")
        f.write(str(q))
        f.write(f"--- {len(N)} points removed --- ")
        f.write(s+"\n")
        f.close()

        M.write(f"Model_{sbox_name}.lp")


    # Dispose
    M.dispose()

    return all_ineqs


def preprocess(all_ineqs, impossible_transitions, n):
    imp_trans_dict = {v: list() for v in impossible_transitions}
    for v in impossible_transitions:
        for i in range(len(all_ineqs)):
            bin_v = format(v, f'0{n}b')
            val = sum(all_ineqs[i][j] * int(bin_v[j]) for j in range(n)) + all_ineqs[i][n]
            if val<0:
                imp_trans_dict[v].append(i)
    return imp_trans_dict


def pick_best_ineqs(P, D, n):

    P_list = list(P)
    m = int(len(P))

    # Create Model
    M = gp.Model()

    # Variables
    d_vars = M.addVars(m, vtype=GRB.BINARY, name="d")

    # Constraints
    for v in D:
        curr_i = set()
        for i in range(m):
            bin_v = format(v, f'0{n}b')
            f_v = sum(P_list[i][j] * int(bin_v[j]) for j in range(n)) + P_list[i][n]
            if f_v < 0:
                curr_i.add(i)
        M.addConstr(gp.quicksum(d_vars[i] for i in curr_i) >= 1)

    # Objective
    M.setObjective(gp.quicksum(d_vars[i] for i in range(m)), GRB.MINIMIZE)

    # Optimize
    M.optimize()

    # Final Inequalities
    final_ineqs_list = list()
    if M.status == GRB.OPTIMAL:
        for i in range(m):
            if int(d_vars[i].x) == 1:
                final_ineqs_list.append(P_list[i])
    else:
        print("No solution found")

    # Dispose
    M.dispose()
    
    return final_ineqs_list

if __name__ == "__main__":

    f = open('../SBOXES/4_bit_sboxes.json')
    # Take S boxes from file input
    data = json.load(f)
    f.close()

    for val in data:
        s = "$"
        for ch in val["name"]:
            if ch == '_':
                s += "\\_"
            else:
                s += ch.upper()
        s +="$"
        # Load S box from data
        name, input_bit_size, output_bit_size, s_box = get_sbox(val)

        # Generate DDT
        DDT = gen_DDT(input_bit_size, output_bit_size, s_box)

        # Get the possible and impossible transitions
        possible_transitions, impossible_transitions = get_transitions(input_bit_size, output_bit_size, DDT) 

        # Intialize variables for this method
        n = input_bit_size+output_bit_size

        start_time = time.time()
        # Final Inequalities
        all_ineqs = gen_functions(possible_transitions, impossible_transitions, n, name)
        all_ineqs_list = list(all_ineqs)
        N = int(len(all_ineqs_list))
        imp_trans_dict = preprocess(all_ineqs_list, impossible_transitions, n)

        imp_trans_funs_list = []
        for v in impossible_transitions:
            imp_trans_funs_list.append(imp_trans_dict[v])

        # Final Inequalities
        Final_inequalities = pick_best_ineqs(all_ineqs, impossible_transitions, n)
        end_time = time.time()

        s = s + " & " + str(len(impossible_transitions)) + " & " + str(len(possible_transitions))
        s += " & " + str(len(Final_inequalities))
        s += " & " + str("{:.3f}".format(end_time - start_time))
        s += "\\\\"
        # print(s)
        file1 = open("new_4_bit_greedy_gen_and_red.txt", "a")
        file1.write(s+"\n")
        file1.close()

        file2 = open(f"{name}.txt", "a")
        for q in Final_inequalities:
            file2.write(str(q)+'\n')
        file2.close()
        
