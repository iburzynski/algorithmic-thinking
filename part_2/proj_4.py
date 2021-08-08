"""
Project #4
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
    scoring = {}
    for char_1 in alpha_copy:
        scoring[char_1] = {}
        for char_2 in alpha_copy:
            if char_1 == char_2 and char_2 != '-':
                score = diag_score
            elif char_1 != '-' and char_2 != '-':
                score = off_diag_score
            else:
                score = dash_score

            scoring[char_1][char_2] = score

    return scoring

def compute_alignment_matrix(seq_x, seq_y, scoring, global_flag):
    """
    Computes an alignment matrix...
    """
    len_m = len(seq_x)
    len_n = len(seq_y)
    if global_flag == False:
        floor = 0
    else:
        floor = -float('inf')

    table = [[0 for dummy_col in range(len_n + 1)] for dummy_row in range(len_m + 1)]
    for idx_i in range(1, len_m + 1):
        table[idx_i][0] = max(table[idx_i - 1][0] + scoring[seq_x[idx_i - 1]]['-'], floor)
    for idx_j in range(1, len_n + 1):
        table[0][idx_j] = max(table[0][idx_j - 1] + scoring['-'][seq_y[idx_j - 1]], floor)
    for idx_i in range(1, len_m + 1):
        for idx_j in range(1, len_n + 1):
            scores = [floor]
            scores.append(table[idx_i - 1][idx_j - 1] + scoring[seq_x[idx_i - 1]][seq_y[idx_j - 1]])
            scores.append(table[idx_i - 1][idx_j] + scoring[seq_x[idx_i - 1]]['-'])
            scores.append(table[idx_i][idx_j - 1] + scoring['-'][seq_y[idx_j - 1]])
            table[idx_i][idx_j] = max(scores)

    return table

def compute_alignment(seq_x, seq_y, scoring, table, global_flag):
    """
    Computes an alignment...
    """
    str_x = ''
    str_y = ''

    if global_flag == False:
        score = max(map(max, table))
        for idx_k, row in enumerate(table):
            for idx_l, col in enumerate(row):
                if col == score:
                    idx_i = idx_k
                    idx_j = idx_l
                    break
    else:
        idx_i = len(seq_x)
        idx_j = len(seq_y)
        score = table[idx_i][idx_j]        

    while idx_i != 0 and idx_j != 0:
        if table[idx_i][idx_j] == table[idx_i - 1][idx_j - 1] + scoring[seq_x[idx_i - 1]][seq_y[idx_j - 1]]:
            str_x = seq_x[idx_i - 1] + str_x
            str_y = seq_y[idx_j - 1] + str_y
            idx_i -= 1
            idx_j -= 1
        elif table[idx_i][idx_j] == table[idx_i - 1][idx_j] + scoring[seq_x[idx_i - 1]]['-']:
            str_x = seq_x[idx_i - 1] + str_x
            str_y = '-' + str_y
            idx_i -= 1
        else:
            str_x = '-' + str_x
            str_y = seq_y[idx_j - 1] + str_y
            idx_j -= 1

        if table[idx_i][idx_j] == 0 and global_flag == False:
            return score, str_x, str_y

    while idx_i != 0:
        str_x = seq_x[idx_i - 1] + str_x
        str_y = '-' + str_y
        idx_i -= 1
    while idx_j != 0:
        str_x = '-' + str_x
        str_y = seq_y[idx_j - 1] + str_y
        idx_j -= 1

    return score, str_x, str_y

def compute_global_alignment(seq_x, seq_y, scoring, a_matrix):
    """
    Wrapper function for computing global alignments 
    (per autograder requirement)
    """
    return compute_alignment(seq_x, seq_y, scoring, a_matrix, global_flag=True)

def compute_local_alignment(seq_x, seq_y, scoring, a_matrix):
    """
    Wrapper function for computing local alignments 
    (per autograder requirement)
    """
    return compute_alignment(seq_x, seq_y, scoring, a_matrix, global_flag=False)
    
### TESTING ###    
scoring = build_scoring_matrix(set(['A', 'C', 'T', 'G']), 10, 4, -6)
test_g = compute_global_alignment('AA', 'TAAT', scoring)
test_l = compute_local_alignment('AA', 'TAAT', scoring)
print(test_g)
print(test_l)