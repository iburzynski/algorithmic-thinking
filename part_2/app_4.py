"""
Provide code and solution for Application 4
"""
import random
import requests
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from proj_4 import compute_alignment_matrix, compute_alignment

# URLs for data files
BASE_DIR = "http://storage.googleapis.com"
ALG_DIR = "/codeskulptor-alg/"
PAM50_URL = f"{BASE_DIR}{ALG_DIR}alg_PAM50.txt"
HUMAN_EYELESS_URL = f"{BASE_DIR}{ALG_DIR}alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = f"{BASE_DIR}{ALG_DIR}alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = f"{BASE_DIR}{ALG_DIR}alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = f"{BASE_DIR}/codeskulptor-assets/assets_scrabble_words3.txt"

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
    word_list = words.split('\n')
    print(f"Loaded a dictionary with {len(word_list)} words")

    return word_list

##########################################################
# Functions for Question 1
##########################################################

def load_data():
    human = read_protein(HUMAN_EYELESS_URL)
    fly = read_protein(FRUITFLY_EYELESS_URL)
    scoring = read_scoring_matrix(PAM50_URL)
    
    return human, fly, scoring

def human_fly_align(human, fly, scoring, randomize=False, mute=True):
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

def question_1(human, fly, scoring, mute=False):
    """
    """
    human, fly, scoring = load_data()

    return human_fly_align(human, fly, scoring, mute=False)

##########################################################
# Functions for Question 2
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

##########################################################
# Functions for Question 4
##########################################################

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

def question_4():
    """
    """
    def normalize_dist(raw_dist):
        """
        Normalizes a distribution (input as a dictionary of key value pairs)
        """
        total = sum(raw_dist.values())

        return {key: (score / total) for key, score in raw_dist.items()} 

    def plot_dist(dist, from_file=True):
        """
        Plots the input distribution and saves as a .png file.
        """
        if from_file and Path.is_file('/data/norm_dist.csv'):
            df = pd.read_csv('data/norm_dist.csv')
        else:
            Path('/data/app_4').mkdir(parents=True, exist_ok=True)
            data = sorted([[score, freq] for score, freq in dist.items()])
            df = pd.DataFrame(data, columns=['Score', 'Frequency'])
            df.set_index('Score', inplace=True)
            df.to_csv('data/app_4/norm_dist.csv')
        print(df.head())
        plt.figure()
        df.plot.bar()
        plt.title(f"Null Distribution for Hypothesis Testing Using {len(data)} Trials")
        plt.xlabel('Scores')
        plt.ylabel('Fraction of Trials')
        sns.despine()

        Path('/plots/app_4').mkdir(parents=True, exist_ok=True)
        plt.savefig(f'plots/app_4/norm_dist_{len(data)}.png')

        return None

    human, fly, scoring = load_data()
    raw_dist = generate_null_distribution(human, fly, scoring, 1000)
    norm_dist = normalize_dist(raw_dist)
    plot_dist(norm_dist)

    return None

##########################################################
# Question 1
##########################################################
human, fly, scoring = load_data()

q1_score, alignment_human, alignment_fly = question_1(human, fly, scoring)

##########################################################
# Question 2
##########################################################

h_corr, f_corr = question_2()

##########################################################
# Question 4
##########################################################

question_4()