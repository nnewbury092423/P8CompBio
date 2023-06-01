#!/usr/bin/env python3
ERROR_PROB = "Stationary probabilities must be between 0 and 1"
ERROR_SEQS = "Invalid sequence file"
PROB_KEYS = ['A', 'C', 'G', 'T']
RATE_KEYS = ['CT', 'AT', 'GT', 'AC', 'CG', 'AG']
import numpy as np
import scipy as sp
from estimate_gtr_pair import gtr_params_pair

def gtr_params(tree, seqs):
    '''
    This function estimates GTR model parameters from a tree and sequences
    :param tree: A TreeSwift Tree object representing the phylogenetic tree
    :param seqs: A dictionary where keys are labels corresponding to the labels of the leaves in ``tree`` and values are sequences (strings)
    :return: A dictionary ``gtr_probs`` storing the GTR stationary probabilities, and a dictionary ``gtr_rates`` storing the GTR transition rates
    '''
    gtr_probs = dict() # keys: {'A', 'C', 'G', 'T'}   values: the corresponding GTR stationary probabilities
    gtr_rates = dict() # keys: {'AC', 'AG', 'AT', 'CG', 'CT', 'GT'}   values: the corresponding GTR transition rates
    # TODO Your code here
    gtr_rates = {'CT':0,'AT':0, 'GT':0, 'AC':0, 'CG':0, 'AG':0}
    gtr_probs = {'A':0, 'C':0, 'G':0, 'T':0}
    count = 0
    for node1 in tree.traverse_postorder():
        if node1.is_leaf():
            for node2 in tree.traverse_postorder():
                if node2.is_leaf() and node1 != node2:
                    d = tree.distance_between(node1,node2)
                    temp_probs,temp_rates = gtr_params_pair(seqs[node1.get_label()],seqs[node2.get_label()],d)
                    for key, value in temp_rates.items():


                        # weight longer branches more than shorter branches to discourage outliers
                        # will normalize later
                        gtr_rates[key] += d*value
                    for key, value in temp_probs.items():
                        gtr_probs[key] += value

    #import pdb; pdb.set_trace()
    # normalize
    norm = gtr_rates['GT']
    for key, value in gtr_rates.items():
        gtr_rates[key] = value/norm

    #normalize
    import pdb; pdb.set_trace()
    probsum = sum(list(gtr_probs.values()))
    for key, value in gtr_probs.items():
        gtr_probs[key] = value/probsum



    
    import pdb; pdb.set_trace()
    #gtr_params_pair(r,s,d)
    return gtr_probs,gtr_rates

def read_FASTA(filename):
    '''
    This function reads a FASTA file from a file and returns a dictionary mapping identifiers to sequences
    '''
    stream = open(filename); seqs = {}; name = None; seq = ''
    for line in stream:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq; stream.close()
    return seqs

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your functions, and outputting the results
    '''
    # parse user arguments
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, type=str, default='stdin', help="Input Tree File (Newick)")
    parser.add_argument('-s', '--seqs', required=True, type=str, help="Multiple Sequence Alignment (FASTA)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output GTR Parameters")
    args = parser.parse_args()

    # load input tree
    from treeswift import read_tree_newick
    if args.tree == 'stdin':
        from sys import stdin
        tree = read_tree_newick(stdin)
    else:
        tree = read_tree_newick(args.tree)
    leaves = {l.label for l in tree.traverse_leaves()}

    # load sequences
    try:
        seqs = read_FASTA(args.seqs)
    except:
        raise ValueError(ERROR_SEQS)
    for k in seqs:
        assert k in leaves, "Sequence ID not in tree: %s" % k

    # run student code and output
    gtr_probs, gtr_rates = gtr_params(tree,seqs)
    for k in PROB_KEYS:
        assert k in gtr_probs, "Missing GTR stationary probability: %s" % k
        if gtr_probs[k] < 0 or gtr_probs[k] > 1:
            raise ValueError(ERROR_PROB)
    for k in RATE_KEYS:
        assert k in gtr_rates, "Missing GTR transition rate: %s" % k
    if args.output == 'stdout':
        from sys import stdout as outfile
    else:
        outfile = open(args.output,'w')
    outfile.write('%f %f %f %f\n' % tuple(gtr_probs[k] for k in PROB_KEYS))
    outfile.write('%f %f %f %f %f %f\n' % tuple(gtr_rates[k] for k in RATE_KEYS))
    outfile.close()
