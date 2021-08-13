"""
Code for Application component of Module 4
Author: Ian Burzynski
"""

# Import helper functions, constants and project functions
import helpers
from constants import PAM50_URL, CONSENSUS_PAX_URL
from proj_4 import compute_alignment_matrix, compute_alignment

# General imports
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
from math import factorial

##########################################################
# Question 1
##########################################################

def question_1(mute=False):
    """
    """
    human, fly, scoring = helpers.load_hf_data()

    return helpers.human_fly_align(human, fly, scoring, mute=mute)

# human, fly, scoring = load_data()
q1_score, alignment_human, alignment_fly = question_1()

##########################################################
# Question 2
##########################################################

def question_2():
    """
    Function for answering Question #2. Returns two percentage values 
    """
    # Define constants
    PAX = helpers.read_protein(CONSENSUS_PAX_URL)
    SCOR = helpers.read_scoring_matrix(PAM50_URL)

    def compare_seq(seq):
        """
        Helper function for comparing sequence to the PAX domain.
        Returns a percentage value (string rounded to 3 significant figures) 
        indicating the correspondence of the two sequences. 
        """
        # Remove dashes from sequence
        seq = seq.replace('-', '')
        a_matrix = compute_alignment_matrix(seq, PAX, SCOR, global_flag=True)
        _, p_align, s_align = compute_alignment(seq, PAX, SCOR, a_matrix, True)
        sequences = zip(p_align, s_align)
        matches = 0
        for p_char, s_char in sequences:
            if p_char == s_char:
                matches += 1

        return f"{matches / len(p_align) * 100:.1f}"

    human, fly, scoring = helpers.load_hf_data()
    _score, h_align, f_align = helpers.human_fly_align(human, fly, scoring)

    h_corr = compare_seq(h_align)
    f_corr = compare_seq(f_align)

    print(f"Human PAX Correspondence: {h_corr}%")
    print(f"Fruitfly PAX Correspondence: {f_corr}% \n")

    return h_corr, f_corr

h_corr, f_corr = question_2()

##########################################################
# Question 3
##########################################################

def question_3():
    """
    """
    def binomial_prob(trials, successes, num_outcomes):
        """
        Helper function to compute the binomial probability of a given number of
        successes out of a given number of trials, given a range of outcomes of 
        size num_outcomes.
        """
        # Compute value of binomial coefficient (n choose k)
        perms = factorial(trials)
        b_coef = perms / (factorial(successes) * factorial(trials - successes))

        p_success = 1 / a_len
        p_fail = (num_outcomes - 1) / num_outcomes
        failures = trials - successes

        return b_coef * (p_success ** successes) * (p_fail ** failures)

    *_seqs, scoring = helpers.load_hf_data()
    a_len = len(scoring) - 1
    _score, h_align, f_align = question_1(mute=True)
    seq_len = len(h_align)
    matches = 0
    for idx, _char in enumerate(h_align):
        if h_align[idx] == f_align[idx]:
            matches += 1

    print(f"Alphabet Length: {a_len}")
    print(f"Alignment Length: {seq_len}")
    print(f"Matching Characters: {matches}")
    print(f"Match Percentage: {matches / seq_len * 100:.1f}")
    
    # Compute the probability of matches being >= the alignment's matches
    match_range = range(matches, seq_len + 1)
    p_total = sum([binomial_prob(seq_len, num, a_len) for num in match_range])
    
    print(f"Probability of {matches} or more matches: {p_total:.1g}\n")

    return p_total

question_3()

##########################################################
# Question 4
##########################################################

def question_4(trials=1000):
    """
    """
    def plot_dist(dist):
        """
        Plots the input distribution and saves as a .png file.
        """
        plt.figure()
        plt.bar(list(dist.keys()), dist.values(), color='#006666', width=1, 
                     edgecolor='#004c4c', linewidth=1)
        plt.title(f"Null Distribution for Hypothesis Testing ({trials} trials)")
        plt.xlabel('Scores')
        plt.tick_params(axis='x', which='major', labelsize=7)
        plt.ylabel('Fraction of Trials')
        sns.despine()

        os.makedirs(f'./plots', exist_ok=True)
        plt.savefig(f'plots/norm_dist_{trials}.png')

        return None

    fname = f'norm_dist_{trials}.pkl'
    dist = helpers.load_dist(fname, trials)
    plot_dist(dist)

    return None

question_4()

##########################################################
# Question 5
##########################################################

def question_5(num_trials=1000):
    """
    """
    # Get the alignment score for human vs. fruitfly sequences 
    q1_score, *_aligns = question_1(mute=True)
    # Load the normal distribution
    fname = f'norm_dist_{num_trials}.pkl'
    dist = helpers.load_dist(fname, num_trials)
    # Reconstruct the trial outcomes from the distribution dictionary
    scores = []
    for score, frq in dist.items():
        scores += [score] * int(frq * num_trials)
    # Compute the mean, std and human vs. fruitfly z-score
    dist_mean = np.mean(scores)
    dist_std = np.std(scores)
    z_score = (q1_score - dist_mean) / dist_std

    print(f"Normal Distribution Mean: {dist_mean:.1f}")
    print(f"Normal Distribution Standard Deviation: {dist_std:.1f}")
    print(f"Human-Fruitfly Alignment Z-score: {z_score:.1f}\n")

    return dist_mean, dist_std, z_score

question_5()

##########################################################
# Question 8
##########################################################

def question_8(itrs, *queries):
    """
    Queries must be formatted as two-element tuples in form ("query", distance), 
    with query as a string and distance as an integer indicating the maximum 
    possible edit distance for matching words to the query.
    """
    time = helpers.time_spell_check('firefly', 2, itrs)
    print(f"Avg. time of check_spelling ({itrs} runs): {time:.2f} sec.\n")

    return helpers.spell_check(*queries, mute=False)

question_8(50, ('humble', 1), ('firefly', 2))

##########################################################
# Question 9
##########################################################

def question_9(itrs, *queries):
    """
    Queries must be formatted as two-element tuples in form ("query", distance), 
    with query as a string and distance as an integer indicating the maximum 
    possible edit distance for matching words to the query.
    """
    time = helpers.time_spell_check('firefly', 2, itrs, fast=True)
    print(f"Avg. time of fast_check_spelling ({itrs} runs): {time:.2f} sec.\n")

    return helpers.spell_check(*queries, fast=True, mute=False), time

question_9(50, ('humble', 1), ('firefly', 2))