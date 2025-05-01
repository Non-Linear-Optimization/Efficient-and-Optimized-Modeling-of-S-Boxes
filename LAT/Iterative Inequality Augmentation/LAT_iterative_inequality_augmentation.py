import numpy as np
import json
import time
from sage.all import *
from itertools import combinations
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


def get_transitions(input_bit_size, output_bit_size, LAT):
    # Get possible and impossible transitions from the LAT
    possible_transitions = list()
    impossible_transitions = list()    
    for i in range(2**input_bit_size):
        for j in range(2**output_bit_size):
            x = bin(i)[2:].zfill(input_bit_size)
            y = bin(j)[2:].zfill(output_bit_size)
            point = [int(bit) for bit in x+y]
            if LAT[i][j] != 0:
                possible_transitions.append(point)
            else:
                impossible_transitions.append(point)
    return possible_transitions, impossible_transitions


def gen_new_ineqs(impossible_transitions, possible_transitions, k):
    P = Polyhedron(vertices = possible_transitions)
    convex_hull = list(P.Hrepresentation())

    candidate_ineqs = set()

    # This is the new part added
    for q in convex_hull:
        candidate_ineqs.add(tuple(list(q)))

    inequality_lists = {q: set() for q in convex_hull}
    for q in convex_hull:
        for point in impossible_transitions:
            sum = q[0]
            for i in range(len(point)):
                sum += (q[i+1]*point[i])
            if sum<0:
                inequality_lists[q].add(tuple(point))

    total_subset_set = set()
    # Iterating over every possible transition
    step = 1
    for point in possible_transitions:
        step += 1

        # For each possible transition, taking all the inequalities it satisfies exactly i.e = 0
        currpoint_inequalities = list()
        for q in convex_hull:
            sum = q[0]
            for i in range(len(point)):
                sum += (q[i+1]*point[i])
            if sum==0:
                currpoint_inequalities.append(q)
        if k <= len(currpoint_inequalities):
            subsets = list(combinations(currpoint_inequalities, k))
            total_imp_transitions_removed = set()

            # Iterating over each subset
            for subset in subsets:
                if subset in total_subset_set:
                    continue
                total_subset_set.add(frozenset(subset))

                # Generating the new inequality
                new_inequality = list(subset[0])
                for i in range(1,len(subset)):
                    for j in range(len(new_inequality)):
                        new_inequality[j] += list(subset[i])[j]

                # Finding the impossible transitions removed by this newly generated inequality
                curr_imp_transitions_removed = set()
                for point in impossible_transitions:
                    sum = new_inequality[0]
                    for i in range(len(point)):
                        sum += (new_inequality[i+1]*point[i])
                    if sum<0:
                        curr_imp_transitions_removed.add(tuple(point))

                # Checking if the impossible transitions removed by this are not in the original removed set
                flag=0
                for q in convex_hull:
                    if len(curr_imp_transitions_removed) == 0:
                        flag=1
                        break
                    if curr_imp_transitions_removed.issubset(inequality_lists[q]):
                        flag=1
                        break
                if flag==0:
                    candidate_ineqs.add(tuple(new_inequality))
        else:
            continue
    return candidate_ineqs, len(convex_hull) 


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


def preprocess(final_ineqs_list, impossible_transitions):
    imp_trans_dict = {tuple(point): list() for point in impossible_transitions}
    for point in impossible_transitions:
        for i in range(len(final_ineqs_list)):
            sum = final_ineqs_list[i][0]
            for j in range(len(point)):
                sum += (final_ineqs_list[i][j+1]*point[j])
            if sum<0:
                imp_trans_dict[tuple(point)].append(i)
    return imp_trans_dict


# In our fina l inequalities, the first term is constant and the terms that follow are coefficients of x1,x2,x3... then y1,y2,y3... respectively
if __name__ == "__main__":

    # All S boxes are loaded into 'data'
    f = open('../SBOXES/5_bit_sboxes.json')
    data = json.load(f)
    f.close()
    # print(data)
    for val in data:
        # Take S box from file input
        name, input_bit_size, output_bit_size, s_box = get_sbox(val)

        # Generate LAT
        LAT = gen_LAT(input_bit_size, output_bit_size, s_box)
        # print(name)
        # print(LAT)
        # print("\n\n")

        # Get possible and impossible transitions
        possible_transitions, impossible_transitions = get_transitions(input_bit_size, output_bit_size, LAT)
        start_time = time.time()
        candidate_ineqs, sage_number = gen_new_ineqs(impossible_transitions, possible_transitions, 2)
        end_time = time.time()
        candidate_ineqs_list = list(candidate_ineqs)
        N = int(len(candidate_ineqs_list))

        imp_trans_dict = preprocess(candidate_ineqs_list, impossible_transitions)
        imp_trans_set = []
        for point in impossible_transitions:
            imp_trans_set.append(imp_trans_dict[tuple(point)])

        # Create a Gurobi model
        model = gp.Model()

        # Create binary variables zi for i = 1 to N
        z = model.addVars(N, vtype=GRB.BINARY, name="z")

        # Objective function: minimize the sum of zi
        model.setObjective(sum(z), GRB.MINIMIZE)

        # Add inequalities to the model
        for eq in imp_trans_set:
            if len(eq)>0:
                model.addConstr(gp.quicksum(z[i] for i in eq) >= 1)
        
        # Optimize the model
        model.setParam('FeasibilityTol', 1e-9)
        model.optimize()
        # model.display()
        # model.printQuality()

        final_ineqs = list()
        if model.status == GRB.OPTIMAL:
            for i in range(N):
                if int(z[i].x)==1:
                    final_ineqs.append(candidate_ineqs_list[i])

        # Dispose of the model
        model.dispose()

        # file = open(f"4-bit_sboxes/{name}_improved_MILP.txt","a")
        # for q in final_ineqs:
        #     file.write(str(list(q))+"\n")
        # file.close()

        f2 = open(f"LAT_IIA_with_time_{input_bit_size}-bit_sboxes.txt","a")
        f2.write("$" + name + "$ & " + str(len(impossible_transitions)) + " & " + str(len(possible_transitions)) + " & " + str(len(final_ineqs)) + " & " + str("{:.1f}".format(end_time - start_time)) + "\n")
        f2.close()
