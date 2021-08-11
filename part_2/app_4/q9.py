import string
import requests

ALPHABET = set(list(string.ascii_lowercase))

BASE_DIR = "http://storage.googleapis.com"
WORD_LIST_URL = f"{BASE_DIR}/codeskulptor-assets/assets_scrabble_words3.txt"

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
    print(f"Loaded a dictionary with {len(word_list)} words")

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

def make_variants(word, dist):
    variants = set([word])
    operations = dist
    while operations > 0:
        vars_copy = variants.copy()
        for word in vars_copy:
            variants.update(insert_letter(word))
            variants.update(delete_letter(word))
            variants.update(substitute_letter(word))
        operations -= 1
    
    word_list = read_words(WORD_LIST_URL)

    return variants.intersection(word_list)

test1 = make_variants('humble', 1)
print(test1)
print(len(test1))