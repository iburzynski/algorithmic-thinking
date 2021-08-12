"""
Helper functions for Application component of Module 4
"""

# Import project code
import proj_4
from constants import HUMAN_EYELESS_URL, FRUITFLY_EYELESS_URL, PAM50_URL, ALPHABET, WORD_LIST_URL

# General imports
import os
import pickle
import random
import requests
import timeit

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

def load_hf_data():
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

    a_matrix = proj_4.compute_alignment_matrix(human, fly_copy, scoring, False)
    align = proj_4.compute_alignment(human, fly_copy, scoring, a_matrix, False)

    if not mute:
        score, h_align, f_align = align
        print(f"Alignment Score: {score}")
        print(f"Human Alignment: {h_align}")
        print(f"Fruitfly Alignment: {f_align}\n")

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
    Normalizes a distribution (input as a dictionary of key value pairs).
    """
    total = sum(raw_dist.values())

    return {key: (score / total) for key, score in raw_dist.items()} 

def save_dist(dist, fname):
    """
    Pickles a distribution dictionary and saves it to the cache folder.
    """
    with open(f'cache/{fname}', 'wb') as file:
        pickle.dump(dist, file, pickle.HIGHEST_PROTOCOL)

    return None

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
        human, fly, scoring = load_hf_data()
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

def check_spelling(query, dist, word_list):
    """
    Inefficient spell checking function that uses global alignment scoring to 
    determine valid word matches. The query word is checked against every word 
    in a provided word list for edit distance less than or equal to the chosen
    value.
    """
    matches = set()
    # Create scoring matrix with diag. score = 2, off-diag. = 1, dash score = 0
    scores = proj_4.build_scoring_matrix(ALPHABET, 2, 1, 0)
    for word in word_list:
        # Get the global alignment score for the query sequence and this word
        matrix = proj_4.compute_alignment_matrix(query, word, scores, True)
        score, *_a = proj_4.compute_alignment(query, word, scores, matrix, True)
        # Check for valid match using formula |Q| + |W| - Score = Edit Distance
        if len(query) + len(word) - score <= dist:
            matches.add(word)

    return matches

def fast_check_spelling(checked_word, dist, word_list):
    """
    Fast implementation of check_spelling using set operations instead of 
    computing alignments.
    """
    # Helper functions for edit operations
    def insert_letter(word):
        """
        Generates a set of all possible insertion variants by inserting each 
        letter of the alphabet at each possible position in the string.
        """
        words = set()
        # Generate all insertion variants for positions up to len(word) - 1
        for pos, _ in enumerate(word):
            for letter in ALPHABET:
                new_word = word[:pos] + letter + word[pos:]
                words.add(new_word)
        # Generate all variants with the letter inserted at the end of the word
        for letter in ALPHABET:
            new_word = word + letter
            words.add(new_word)

        return words

    def delete_letter(word):
        """
        Generate a set of all possible deletion variants by deleting each letter 
        of the string.
        """
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