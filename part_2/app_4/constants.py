"""
Constants for Application component of Module 4
"""

import string

# Constants for file paths
BASE_DIR = "http://storage.googleapis.com"
ALG_DIR = "/codeskulptor-alg/"
PAM50_URL = f"{BASE_DIR}{ALG_DIR}alg_PAM50.txt"
HUMAN_EYELESS_URL = f"{BASE_DIR}{ALG_DIR}alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = f"{BASE_DIR}{ALG_DIR}alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = f"{BASE_DIR}{ALG_DIR}alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = f"{BASE_DIR}/codeskulptor-assets/assets_scrabble_words3.txt"
DATA_PATH = "data/app_4/"
PLOTS_PATH = "plots/app_4/"

# Alphabet set for spell check questions
ALPHABET = set(list(string.ascii_lowercase))