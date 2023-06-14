[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/VzwJPzPo)
# ECE 208 Homework 8: Estimating GTR Parameters
In this homework assignment, we will continue studying models of DNA evolution, namely the [General Time Reversible (GTR) model](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#GTR_model_(Tavar%C3%A9_1986)).<sup>1</sup> In the previous assignment, we tried to compute the (log-)likelihood of observing DNA sequences *S* given the GTR model parameters and given a tree *T* with known branch lengths and topology. Now, we will try to infer the maximum-likelihood GTR model parameters given DNA sequences *S* and tree *T*.

## Estimating GTR Parameters from a Pair of Sequences
We will start with an easier sub-problem: given two sequences *r* and *s* and the distance between them (in unit of expected number of per-site mutations), what are the maximum-likelihood GTR model parameters?

### Input
[estimate_gtr_pair.py](estimate_gtr_pair.py) will take as input the following items:
* `-s`: The path to a [FASTA file](https://en.wikipedia.org/wiki/FASTA_format) containing two sequences *r* and *s*
* `-d`: The pairwise distance between *r* and *s*
* `-o`: The output file, which will contain the the GTR parameters (see below) (default: standard output)

You will only be modifying the `gtr_params_pair(r, s, d)` function, which takes as input the following parameters:
* `r`: A sequence of length *k*
* `s`: A sequence of length *k*
* `d`: The pairwise distance between *r* and *s*

### Output
[estimate_gtr_pair.py](estimate_gtr_pair.py) will output the estimated GTR model parameters in the following format (an example can be found in [example/gtr_params.txt](example/gtr_params.txt)):

```
prob_A prob_C prob_G prob_T
rate_CT rate_AT rate_GT rate_AC rate_CG rate_AG
```

The `gtr_params_pair(r, s, d)` function will return two dictionaries, `gtr_probs` and `gtr_rates`, containing the following items:

```python
# Stationary Frequencies
gtr_probs['A']  # stationary frequency of A
gtr_probs['C']  # stationary frequency of C
gtr_probs['G']  # stationary frequency of G
gtr_probs['T']  # stationary frequency of T

# Transition Rates
gtr_rates['AC'] # transition rate of A <-> C
gtr_rates['AG'] # transition rate of A <-> G
gtr_rates['AT'] # transition rate of A <-> T
gtr_rates['CG'] # transition rate of C <-> G
gtr_rates['CT'] # transition rate of C <-> T
gtr_rates['GT'] # transition rate of G <-> T
```

## Estimating GTR Parameters from Multiple Sequences
Now, we will generalize the previous step to try to estimate GTR model parameters given multiple sequences *S* and a tree *T* with known branch lengths and topology. **Note that there is no single solution, and you need to come up with your own *heuristic* solution!** 

### Input
[estimate_gtr.py](estimate_gtr.py) will take as input the following items:
* `-t`: A phylogenetic tree *T* (either binary or multifurcating) in the [Newick format](https://en.wikipedia.org/wiki/Newick_format)
* `-s`: The path to a [FASTA file](https://en.wikipedia.org/wiki/FASTA_format) containing a multiple sequence alignment *S* such that each sequence in *S* corresponds to a leaf in *T*
* `-o`: The output file, which will contain the the GTR parameters (same format as above) (default: standard output)

You will only be modifying the `gtr_params(tree, seqs)` function, which takes as input the following parameters:
* `tree`: A rooted [TreeSwift](https://github.com/niemasd/TreeSwift) `Tree` object (either binary or multifurcating)
* `seqs`: A dictionary where keys are labels corresponding to the labels of the leaves in `tree` and values are sequences (strings)

### Output
[estimate_gtr.py](estimate_gtr.py) will output the estimated GTR model parameters in the same format as described above. The `gtr_params(tree, seqs)` function will return two dictionaries, `gtr_probs` and `gtr_rates`, in the same format as described above.

### Examples
We have provided example datasets of sequence pairs ([example/pair](example/pair)) and multiple sequences ([example/multi](example/multi)). The input sequence datasets in the FASTA format are `*.fas`. For pairs of sequences, their corresponding distances are `*.dist.txt`. For multiple sequences, they all correspond to the tree [example/tree.nwk](example/tree.nwk). For all datasets, the true GTR model parameters are in [example/gtr_params.txt](example/gtr_params.txt), and the GTR parameters we inferred from them are `*.est_gtr.txt` (in the same format as described above). Short sequences (*k* = 200) will be in the `short` folders, and long sequences (*k* = 2,000) will be in the `long` folders.

## Homework Deliverables
* **[estimate_gtr_pair.py](estimate_gtr_pair.py):** Your Python 3 code for estimating GTR parameters from a pair of sequences and a distance
    * Usage: `python3 estimate_gtr_pair.py -s <seqs> -d <distance> -o <out_file>`
    * We have provided starter code in the file, but you must fill out the `gtr_params_pair(r, s, d)` function
        * **Note:** Feel free to add import statements to NumPy and SciPy within this function as you need!
* **[estimate_gtr.py](estimate_gtr.py):** Your Python 3 code for estimating GTR parameters from multiple sequences and a tree
    * Usage: `python3 estimate_gtr.py -t <tree> -s <seqs> -o <out_file>`
    * We have provided starter code in the file, but you must fill out the `gtr_params(tree, seqs)` function
        * **Note:** Feel free to add import statements to NumPy and SciPy within this function as you need!
* **[writeup.pdf](writeup.pdf):** A write-up that includes the following pieces of information:
    1. Your derivation of the equations for estimating the maximum-likelihood GTR parameters from a pair of sequences and a distance
    2. A description of your algorithm to estimate GTR parameters from multiple sequences and a tree

## Grade Breakdown (100 Points)
* **[writeup.pdf](writeup.pdf): 20 Points**
    * *Equations for Estimating from a Pair of Sequences:* 10 points
    * *Algorithm for Estimating from Multiple Sequences:* 10 points
* **[estimate_gtr_pair.py](estimate_gtr_pair.py): 40 Points**
    * *Accuracy:* 40 points
        * 5 replicates of *k* = 200, each worth 4 points
        * 5 replicates of *k* = 2,000, each worth 4 points
* **[estimate_gtr.py](estimate_gtr.py): 40 Points**
    * *Accuracy:* 40 points
        * 5 replicates of *k* = 200, each worth 4 points
        * 5 replicates of *k* = 2,000, each worth 4 points

## References
1. Tavaré S. Some probabilistic and statistical problems in the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences*, 17:57–86, 1986.
