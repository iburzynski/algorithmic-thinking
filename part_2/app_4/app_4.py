"""
Provide code and solution for Application 4
"""
import random
import requests
import os
import pickle
import string
import timeit
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from proj_4 import build_scoring_matrix, compute_alignment_matrix, compute_alignment

# URLs for data files
BASE_DIR = "http://storage.googleapis.com"
ALG_DIR = "/codeskulptor-alg/"
PAM50_URL = f"{BASE_DIR}{ALG_DIR}alg_PAM50.txt"
HUMAN_EYELESS_URL = f"{BASE_DIR}{ALG_DIR}alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = f"{BASE_DIR}{ALG_DIR}alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = f"{BASE_DIR}{ALG_DIR}alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = f"{BASE_DIR}/codeskulptor-assets/assets_scrabble_words3.txt"

DATA_PATH = "data/app_4/"
PLOTS_PATH = "plots/app_4/"
ALPHABET = set(list(string.ascii_lowercase))

###############################################
# Provided Functions
###############################################

def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.  

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = requests.get(filename)
    scoring_lines = scoring_file.text.split('\n')
    ykeychars = scoring_lines[0].split()

    for line in scoring_lines[1:]:
        vals = line.split()
        xkey = vals.pop(0)
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            scoring_dict[xkey][ykey] = int(val)

    return scoring_dict

def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    protein_file = requests.get(filename)
    protein_seq = protein_file.text
    protein_seq = protein_seq.rstrip()

    return protein_seq

def load_data():
    """
    Loads and returns the human and fruitfly proteins and the PAM50 scoring 
    matrix.
    """
    human = read_protein(HUMAN_EYELESS_URL)
    fly = read_protein(FRUITFLY_EYELESS_URL)
    scoring = read_scoring_matrix(PAM50_URL)
    
    return human, fly, scoring

def human_fly_align(human, fly, scoring, randomize=False, mute=True):
    """
    Computes the local alignment between human and fruitfly sequences.
    Fruitfly sequence can be randomized.
    Mute flag turns off print statements.
    """
    fly_copy = fly
    # Randomization code for Question #4
    if randomize:
        fly_split = list(fly_copy)
        random.shuffle(fly_split)
        fly_copy = ''.join(fly_split)

    a_matrix = compute_alignment_matrix(human, fly_copy, scoring, global_flag=False)
    align = compute_alignment(human, fly_copy, scoring, a_matrix, global_flag=False)

    if not mute:
        score, h_align, f_align = align
        print(f"Alignment Score: {score}")
        print(f"Human Alignment: {h_align}")
        print(f"Fruitfly Alignment: {f_align}")

    return align

def generate_null_distribution(seq_x, seq_y, scoring, num_trials):
    """
    """
    scoring_distribution = {}
    
    for _trial in range(num_trials):
        score, *_seqs = human_fly_align(seq_x, seq_y, scoring, randomize=True)

        if score in scoring_distribution:
            scoring_distribution[score] += 1
        else:
            scoring_distribution[score] = 1

    return scoring_distribution

def normalize_dist(raw_dist):
    """
    Normalizes a distribution (input as a dictionary of key value pairs)
    """
    total = sum(raw_dist.values())

    return {key: (score / total) for key, score in raw_dist.items()} 

def save_dist(dist, fname):
    """
    Pickles a distribution dictionary and saves it to the cache folder.
    """
    with open(f'cache/{fname}', 'wb') as file:
        pickle.dump(dist, file, pickle.HIGHEST_PROTOCOL)

def load_dist(fname, num_trials=1000):
    """
    Loads a pickled distribution dictionary from cache if it exists. Otherwise, 
    generates a new null distribution for the specified number of trials and 
    saves it to cache as a pickle file.
    """
    if os.path.exists(f'./cache/{fname}'):
        with open(f'cache/{fname}', 'rb') as file:
            return pickle.load(file)
    else:
        # Load data and generate null distribution
        human, fly, scoring = load_data()
        raw_dist = generate_null_distribution(human, fly, scoring, num_trials)
        dist = normalize_dist(raw_dist)

        # Pickle distribution for faster output in future runs
        os.makedirs(f'./cache', exist_ok=True)
        save_dist(dist, fname)

        return dist

def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = requests.get(filename)
    
    # read in files as string
    words = word_file.text
    
    # template lines and solution lines list of line string
    word_list = set(words.split('\n'))
    # print(f"Loaded a dictionary with {len(word_list)} words")

    return word_list

def insert_letter(word):
    words = set()
    for pos, _ in enumerate(word):
        for letter in ALPHABET:
            new_word = word[:pos] + letter + word[pos:]
            words.add(new_word)
    for letter in ALPHABET:
        new_word = word + letter
        words.add(new_word)

    return words

def delete_letter(word):
    words = set()
    for pos, _ in enumerate(word):
        new_word = word[:pos] + word[pos + 1:]
        words.add(new_word)

    return words

def substitute_letter(word):
    words = set()
    for pos, _ in enumerate(word):
        for letter in ALPHABET:
            new_word = word[:pos] + letter + word[pos + 1:]
            words.add(new_word)

    return words

def check_spelling(checked_word, dist, word_list):
    """
    """
    matches = set()
    scoring = build_scoring_matrix(ALPHABET, 2, 1, 0)
    for word in word_list:
        a_matrix = compute_alignment_matrix(checked_word, word, scoring, True)
        score, *_alignments = compute_alignment(checked_word, word, scoring, a_matrix, True)
        if len(checked_word) + len(word) - score <= dist:
            matches.add(word)

    return matches

def fast_check_spelling(checked_word, dist, word_list):
    """
    Fast implementation of check_spelling using set operations instead of 
    computing alignments.
    """
    variants = set([checked_word])
    operations = dist
    while operations > 0:
        vars_copy = variants.copy()
        for word in vars_copy:
            variants.update(insert_letter(word))
            variants.update(delete_letter(word))
            variants.update(substitute_letter(word))
        operations -= 1

    return variants.intersection(word_list)

def spell_check(*queries, fast=False, mute=True):
    if fast:
        func = fast_check_spelling
    else:
        func = check_spelling

    word_list = read_words(WORD_LIST_URL)
    results = {}
    for query, dist in queries:
        matches = func(query, dist, word_list)
        if not mute:
            print(f"Words within Edit Distance {dist} of '{query}':")
            print(matches, '\n')
        results[query] = matches

    return results

def time_spell_check(query, dist, iterations, fast=False):
    """
    """
    if fast:
        func = fast_check_spelling
    else:
        func = check_spelling

    word_list = read_words(WORD_LIST_URL)    

    def closure():
        """
        Closure for running timeit test with provided arguments.
        """
        return func(query, dist, word_list)

    return timeit.timeit(closure, number = iterations) / iterations

##########################################################
# Question 1
##########################################################

def question_1():
    """
    """
    human, fly, scoring = load_data()

    return human_fly_align(human, fly, scoring, mute=False)

# human, fly, scoring = load_data()
q1_score, alignment_human, alignment_fly = question_1()

##########################################################
# Question 2
##########################################################

def question_2():
    """
    Function for answering Question #2. Returns two percentage values 
    """
    # define constants
    PAX = read_protein(CONSENSUS_PAX_URL)
    SCOR = read_scoring_matrix(PAM50_URL)
    def compare_seq(seq):
        """
        Helper function for comparing sequence to the PAX domain.
        Returns a percentage value (string rounded to 3 significant figures) 
        indicating the correspondence of the two sequences. 
        """
        # remove dashes from sequence
        seq = seq.replace('-', '')
        a_matrix = compute_alignment_matrix(seq, PAX, SCOR, global_flag=True)
        _, p_align, s_align = compute_alignment(seq, PAX, SCOR, a_matrix, True)
        sequences = zip(p_align, s_align)
        matches = 0
        for p_char, s_char in sequences:
            if p_char == s_char:
                matches += 1

        return f"{matches / len(p_align) * 100:.1f}"

    human, fly, scoring = load_data()
    _, h_align, f_align = human_fly_align(human, fly, scoring)

    h_corr = compare_seq(h_align)
    f_corr = compare_seq(f_align)

    print(f"Human PAX Correspondence: {h_corr}%")
    print(f"Fruitfly PAX Correspondence: {f_corr}%")

    return h_corr, f_corr

h_corr, f_corr = question_2()

##########################################################
# Question 4
##########################################################

def question_4(num_trials=1000):
    """
    """

    def plot_dist(dist):
        """
        Plots the input distribution and saves as a .png file.
        """
        plt.figure()
        plt.bar(list(dist.keys()), dist.values(), color='#006666', width=1, edgecolor='#004c4c', linewidth=1)
        plt.title(f"Null Distribution for Hypothesis Testing Using {num_trials} Trials")
        plt.xlabel('Scores')
        plt.tick_params(axis='x', which='major', labelsize=7)
        plt.ylabel('Fraction of Trials')
        sns.despine()
        # Save the plot
        os.makedirs(f'./plots', exist_ok=True)
        plt.savefig(f'plots/norm_dist_{num_trials}.png')

        return None

    fname = f'norm_dist_{num_trials}.pkl'
    dist = load_dist(fname, num_trials)
    plot_dist(dist)

    return None

question_4()

##########################################################
# Question 5
##########################################################

def question_5(num_trials=1000):
    # Get the alignment score for human vs. fruitfly sequences 
    q1_score, *_alignments = question_1()
    # Load the normal distribution
    fname = f'norm_dist_{num_trials}.pkl'
    dist = load_dist(fname, num_trials)
    scores = [([score] * int(freq * num_trials)) for score, freq in dist.items()]
    flat_scores = [item for sublist in scores for item in sublist]
    # Compute the mean, std and human vs. fruitfly z-score
    dist_mean = np.mean(flat_scores)
    dist_std = np.std(flat_scores)
    z_score = (q1_score - dist_mean) / dist_std

    print(dist_mean, dist_std, z_score)

question_5()

##########################################################
# Question 8
##########################################################

def question_8(iterations, *queries):
    """
    Queries must be formatted as two-element tuples in form ("query", distance), 
    with query as a string and distance as an integer indicating the maximum 
    possible edit distance for matching words to the query.
    """
    time = time_spell_check('firefly', 2, iterations)
    print(f"Completed {iterations * len(queries)} tests of check_spelling.")
    print(f"Average running time: {time} seconds.")

    return spell_check(*queries, mute=False)

question_8(50, ('humble', 1), ('firefly', 2))

##########################################################
# Question 9
##########################################################

def question_9(iterations, *queries):
    """
    Queries must be formatted as two-element tuples in form ("query", distance), 
    with query as a string and distance as an integer indicating the maximum 
    possible edit distance for matching words to the query.
    """
    time = time_spell_check('firefly', 2, iterations, fast=True)
    print(f"Completed {iterations} tests of fast_check_spelling.")
    print(f"Average running time: {time} seconds.")

    return spell_check(*queries, fast=True, mute=False), time

question_9(50, ('humble', 1), ('firefly', 2))