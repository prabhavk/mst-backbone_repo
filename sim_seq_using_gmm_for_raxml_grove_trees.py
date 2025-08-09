from ete3 import Tree
import numpy as np
from numpy.random import uniform, choice, dirichlet
import sys
_ = sys.argv[0]
DNA=["A","C","G","T"]


def SampleFromDistribution(rootProb,length):    
    seq=""    
    for _ in range(length):
        u = uniform()
        cumSum=0
        for base in range(4):
            cumSum += rootProb[base]
            if cumSum > u:
                break
        seq+=DNA[base]
    return seq

def GenerateEvolveCharFunction(P):
    
    def EvolveChar(char):
        row = DNA.index(char)
        P_row = P[row,:]
        cumSum = 0
        u = uniform()
        for col in range(4):
            cumSum += P_row[col]
            if cumSum > u:
                break
        return DNA[col]
    
    return EvolveChar


def simulate_probs_gc_range(alpha=(1,1,1,1), gc_min=0.13, gc_max=0.75, n=10):
    """
    Simulates nucleotide probabilities (A,C,G,T) from a Dirichlet distribution,
    filtering to ensure pG+pC is within [gc_min, gc_max].

    alpha: Dirichlet concentration parameters (default = uniform Dirichlet(1,1,1,1))
    gc_min, gc_max: acceptable GC content bounds
    n: number of samples to return
    """
    results = []
    while len(results) < n:
        p = np.random.dirichlet(alpha)
        pA, pC, pG, pT = p
        gc_content = pG + pC
        if gc_min <= gc_content <= gc_max:
            results.append((pA, pC, pG, pT, gc_content))
    return results

def joint_and_transition(Px, t, alpha_off=1.0):
    """
    Args:
        Px : array-like of shape (4,), marginal distribution for parent state x.
        t  : branch length in [0, 1]; diagonal mass = 1 - t.
        alpha_off : Dirichlet concentration for off-diagonal allocations.
        
    Returns:
        Pxy : (4,4) joint probability table, sum=1, trace=1-t.
        Py_given_x : (4,4) transition probability table, rows sum to 1.
    """
    Px = np.asarray(Px, dtype=float)
    if not np.isclose(Px.sum(), 1.0):
        raise ValueError("Px must sum to 1.")
    if np.any(Px <= 0):
        raise ValueError("Px must have positive entries (for division).")
    if not (0.0 <= t <= 1.0):
        raise ValueError("t must be between 0 and 1.")

    # Step 1: allocate joint probabilities
    Pxy = np.zeros((4, 4), dtype=float)

    # Diagonals: (1 - t) mass, proportional to Px
    if 1 - t > 0:
        np.fill_diagonal(Pxy, (1 - t) * Px)

    # Off-diagonals: t mass, split within each row via Dirichlet
    if t > 0:
        for i in range(4):
            off_mass = t * Px[i]  # total off-diagonal mass for row i
            weights = np.random.dirichlet([alpha_off] * 3)
            off_idx = [j for j in range(4) if j != i]
            for w, j in zip(weights, off_idx):
                Pxy[i, j] = off_mass * w

    # Sanity check
    assert np.isclose(Pxy.sum(), 1.0), "Joint probabilities must sum to 1"
    assert np.isclose(np.trace(Pxy), 1.0 - t), "Trace must be 1 - t"

    # Step 2: compute conditional probabilities P(y|x) = P(x,y) / P(x)
    Py_given_x = Pxy / Px[:, None]

    # Sanity: rows should sum to 1
    assert np.allclose(Py_given_x.sum(axis=1), 1.0), "Rows of P(y|x) must sum to 1"

    return Pxy, Py_given_x

# P_j, P_t = joint_and_transition(pi["Root"],0.09270729298214177)
# pi_child = P_j.sum(axis=1)

# gc_min = 13
# gc_max = 75

def simulate_sequences_for_raxml_grove_tree(tree_id):
    selected_tree_file_path = "/home/pk/projects/RAxMLGrove/selected_trees/" + tree_id + "/tree_best.newick"
    t = Tree(selected_tree_file_path, format=1)
    gmm_file_path = "/home/pk/projects/RAxMLGrove/selected_trees/" + tree_id + "/simulation_parameters.gmm"
    seq_len = 10000
    
    sequences_file_path = "/home/pk/projects/RAxMLGrove/selected_trees/" + tree_id + "/sequence_len_"+str(seq_len)+".fasta"
    sequences_file = open(sequences_file_path,"w")
    samples = simulate_probs_gc_range(alpha=(1,1,1,1), n=1)
    gmm_file = open(gmm_file_path,"w")
    gmm_file.write("transition matrix for edge from parent_name to child_name\n")
    gmm_file.write("p(A|A) p(C|A) p(G|A) p(T|A) p(A|C) p(C|C) p(G|C) p(T|C) p(A|G) p(C|G) p(G|G) p(T|G) p(A|T) p(C|T) p(G|T) p(T|T)\n")
    pi = {}
    for pA, pC, pG, pT, gc in samples:
        pi["root"] = [pA, pC, pG, pT]
    sequences = {}
    transition_matrices = {}
    sequences["root"] = SampleFromDistribution(pi["root"], seq_len)
    t = Tree(selected_tree_file_path, format=1)
        
    # h_ind = 0
    int_id = 1
    for node in t.traverse("preorder"):    
        pi["root"] = [pA, pC, pG, pT]
        if node.is_root():
            if not node.name:
                node.name = "root"
        else:
            if not node.name:
                node.name = f"h_{int_id}"
                int_id += 1
            parent_name = node.up.name
            gmm_file.write("transition matrix for edge from " + node.up.name +" to " + node.name+"\n")
            branch_length = node.dist
            if (branch_length < 0 or branch_length > 1):
                print("weird branch_length")
                print(branch_length)
                print("----------------------")
                exit
            pi_p = pi[parent_name]
            P_j, P_t = joint_and_transition(pi_p,branch_length)
            transition_matrices[node.name] = P_t
            result_str = ""
            for i in range(4):
                for j in range(4):
                    result_str += f"{round(P_t[i][j],ndigits=6)} "
            
            result_str = result_str.strip()
            gmm_file.write(result_str+ "\n")
            pi[node.name] = P_j.sum(axis=1)
            EvolveChar = GenerateEvolveCharFunction(P_t)
            parent_seq = sequences[parent_name]
            sequences[node.name] = "".join(map(EvolveChar,parent_seq))
            if (not node.name.startswith("h_")):
                sequences_file.write(">" + node.name + "\n" + sequences[node.name]+"\n")

    gmm_file.write("probability at root root is p(A) p(C) p(G) p(T)\n")
    result_str = ""
    for i in range(4):
        result_str += f"{pi["root"][i]} "
    result_str = result_str.strip()
    gmm_file.write(result_str+"\n")
    
# tree_id = sys.argv[1]
# selected_tree_file_path = "/home/pk/projects/RAxMLGrove/trees/" + tree_id + "/tree_best.newick"

# t = Tree(selected_tree_file_path, format=1) 
import subprocess as sub
path_to_grove = "/home/pk/projects/RAxMLGrove/"
selected_trees_file_path = "/home/pk/projects/RAxMLGrove/selected_grove_tree_ids_scalability"

dna_trees_file_path = "/home/pk/projects/RAxMLGrove/all_dna_trees"

dna_trees_file = open(dna_trees_file_path,"r")
dna_trees_list = []
for line in dna_trees_file:
    dna_trees_list.append(line)

dna_trees_file.close()

print(len(dna_trees_list))
selected_trees_file = open(selected_trees_file_path,"r")
# tree_ids = []
line_no = 0
for line in selected_trees_file:    
    # tree_ids.append(line.split(",")[0])    
    break
    tree_id = line.split(",")[0]
    print(line_no, tree_id)
    sub.call("cp -r " + path_to_grove + "trees/" + tree_id + " " + path_to_grove + "selected_trees/", shell=True)
    simulate_sequences_for_raxml_grove_tree(tree_id)
    line_no += 1
    
    

# print(tree_ids)
selected_trees_file.close()


