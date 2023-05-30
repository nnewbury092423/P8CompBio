#!/usr/bin/env python3
ERROR_PROB = "Stationary probabilities must be between 0 and 1"
ERROR_SEQS = "Invalid sequence file"
ERROR_TWO_SEQS = "Sequence file must have exactly 2 sequences"
PROB_KEYS = ['A', 'C', 'G', 'T']
RATE_KEYS = ['CT', 'AT', 'GT', 'AC', 'CG', 'AG']

def gtr_params_pair(r, s, d):
    '''
    This function estimates GTR model parameters from a tree and sequences
    :param r: A sequence of length ``k``
    :param s: A sequence of length ``k``
    :param d: The pairwise distance between ``r`` and ``s``
    :return: A dictionary ``gtr_probs`` storing the GTR stationary probabilities, and a dictionary ``gtr_rates`` storing the GTR transition rates
    '''
    gtr_probs = dict() # keys: {'A', 'C', 'G', 'T'}   values: the corresponding GTR stationary probabilities
    gtr_rates = dict() # keys: {'AC', 'AG', 'AT', 'CG', 'CT', 'GT'}   values: the corresponding GTR transition rates
    # TODO Your code here
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
    parser.add_argument('-s', '--seqs', required=True, type=str, help="Sequences (FASTA)")
    parser.add_argument('-d', '--distance', required=True, type=str, help="Pairwise Distance")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output GTR Parameters")
    args = parser.parse_args()

    # load sequences
    try:
        seqs = read_FASTA(args.seqs)
    except:
        raise ValueError(ERROR_SEQS)
    if len(seqs) != 2:
        raise ValueError(ERROR_TWO_SEQS)
    r,s = seqs.values()

    # load distance
    from os.path import isfile
    if isfile(args.distance):
        f = open(args.distance); d = float(f.read()); f.close()
    else:
        d = float(args.distance)

    # run student code and output
    gtr_probs, gtr_rates = gtr_params_pair(r,s,d)
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
