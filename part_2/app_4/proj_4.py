"""
Functions for Project component of Module 4
Author: Ian Burzynski
"""

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Prepares a nested dictionary from a set of letters (alphabet) and 
    three scores (diag_score: score for pair of matching letters, 
    off_diag_score: score for pair of non-matching letters, dash_score: score 
    for a pair that includes a dash).
    The ouput dictionary has items of the form {char_1: {char_2a: score_a, 
    char_2b: score_b, ..., char_2n: score_n}}, with each score representing 
    the score for the pair of char_1 and char_2.
    """
    alpha_copy = alphabet.copy()
    alpha_copy.add('-')
    # Create the outer dictionary
    scoring = {}
    for char_1 in alpha_copy:
        # Add the character to outer dictionary with an empty inner dictionary
        scoring[char_1] = {}
        # Populate inner dictionary with scores for each possible character pair
        for char_2 in alpha_copy:
            # Choose diagonal score if characters are the same and are letters
            if char_1 == char_2 and char_2 != '-':
                score = diag_score
            # Choose off-diagonal score if characters are different letters
            elif char_1 != '-' and char_2 != '-':
                score = off_diag_score
            # Choose dash score if either character is a dash 
            else:
                score = dash_score
            # Add the chosen score to the inner dictionary
            scoring[char_1][char_2] = score

    return scoring

def compute_alignment_matrix(seq_x, seq_y, scoring, global_flag):
    """
    Computes an alignment matrix for two string sequences using dynamic 
    programming.
    """
    len_x = len(seq_x)
    len_y = len(seq_y)
    # Set the minimum score depending on which type of alignment is computed
    if global_flag == False:
        floor = 0
    else:
        floor = -float('inf')
    # Create a zero matrix based on sequence dimensions
    matrix = [[0 for _col in range(len_y + 1)] for _row in range(len_x + 1)]
    # Populate the dash score cells
    for idx_x in range(1, len_x + 1):
        val = max(matrix[idx_x-1][0] + scoring[seq_x[idx_x - 1]]['-'], floor)
        matrix[idx_x][0] = val
    for idx_y in range(1, len_y + 1):
        val = max(matrix[0][idx_y - 1] + scoring['-'][seq_y[idx_y - 1]], floor)
        matrix[0][idx_y] = val
    # Populate all other cells with the max score of three possible pairings
    for idx_x in range(1, len_x + 1):
        for idx_y in range(1, len_y + 1):
            # Create list of possible score values, including floor value
            scores = [floor]
            # Option 1: xy pair
            x_y = scoring[seq_x[idx_x - 1]][seq_y[idx_y - 1]]
            scores.append(matrix[idx_x - 1][idx_y - 1] + x_y)
            # Option 2: x- pair
            x_dash = scoring[seq_x[idx_x - 1]]['-']
            scores.append(matrix[idx_x - 1][idx_y] + x_dash)
            # Option 3: -y pair
            dash_y = scoring['-'][seq_y[idx_y - 1]]
            scores.append(matrix[idx_x][idx_y - 1] + dash_y)
            # Select the maximum score from the list and assign to cell
            matrix[idx_x][idx_y] = max(scores)

    return matrix

def compute_alignment(seq_x, seq_y, scores, amatrix, global_flag):
    """
    Constructs the optimal global or local alignment for two sequences by 
    performing a traceback over their alignment matrix.
    """
    def find_max_score():
        """
        Helper function for finding the max score and starting indices for 
        computing local alignments.
        """
        # Find the maximum score in the alignment matrix
        score = max(map(max, amatrix))
        # Find the corresponding indices in the matrix
        for idx_x, row in enumerate(amatrix):
            for idx_y, col in enumerate(row):
                if col == score:
                    # Return max score and starting indices when found
                    return score, idx_x, idx_y

    # Initialize empty strings to for reconstructing the optimal alignment
    str_x = ''
    str_y = ''

    # Local alignments: set score and indices from the max score in the matrix
    if global_flag == False:
        score, idx_x, idx_y = find_max_score()
    # Global alignments: set score and indices from the last cell in the matrix
    else:
        idx_x = len(seq_x)
        idx_y = len(seq_y)
        score = amatrix[idx_x][idx_y]

    # Step backwards through sequences using scoring matrix to reconstruct them
    while idx_x != 0 and idx_y != 0:
        # Character pair score values for identifying xy and x- pairs
        x_y =  scores[seq_x[idx_x - 1]][seq_y[idx_y - 1]]
        x_dash = scores[seq_x[idx_x - 1]]['-']
        # If character pair is xy:
        if amatrix[idx_x][idx_y] == amatrix[idx_x - 1][idx_y - 1] + x_y:
            # Prepend x and y characters to the respective sequences
            str_x = seq_x[idx_x - 1] + str_x
            str_y = seq_y[idx_y - 1] + str_y
            # Decrement both indices
            idx_x -= 1
            idx_y -= 1
        # If character pair is x-:
        elif amatrix[idx_x][idx_y] == amatrix[idx_x - 1][idx_y] + x_dash:
            # Prepend x to sequence x and a dash to sequence y
            str_x = seq_x[idx_x - 1] + str_x
            str_y = '-' + str_y
            # Decrement only x-index
            idx_x -= 1
        # If pair is -y:
        else:
            # Prepend a dash to sequence x and y to sequence y
            str_x = '-' + str_x
            str_y = seq_y[idx_y - 1] + str_y
            # Decrement only y-index
            idx_y -= 1
        # For local alignments, return values as soon as a zero value is found
        if amatrix[idx_x][idx_y] == 0 and global_flag == False:

            return score, str_x, str_y

    # Prepend dashes to seq_y for any remaining characters in sequence x
    while idx_x != 0:
        str_x = seq_x[idx_x - 1] + str_x
        str_y = '-' + str_y
        idx_x -= 1
    # Prepend dashes to seq_x for any remaining characters in sequence y
    while idx_y != 0:
        str_x = '-' + str_x
        str_y = seq_y[idx_y - 1] + str_y
        idx_y -= 1

    return score, str_x, str_y

# Functions for testing alignments in autograder (not used in application)
def compute_global_alignment(seq_x, seq_y, scoring, a_matrix):
    """
    Wrapper function for computing global alignments (autograder requirement)
    """
    return compute_alignment(seq_x, seq_y, scoring, a_matrix, global_flag=True)

def compute_local_alignment(seq_x, seq_y, scoring, a_matrix):
    """
    Wrapper function for computing local alignments (autograder requirement)
    """
    return compute_alignment(seq_x, seq_y, scoring, a_matrix, global_flag=False)