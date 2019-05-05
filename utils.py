def calc_first_occ(s):
    """ calculate the first occurance of a letter in sorted string s """
    # s - is the bwt transformed string
    A = {}  # letter count
    for i, c in enumerate(s):
        if A.get(c):
            A[c] += 1
        else:
            A[c] = 1

    # sort the letters
    letters = sorted(A.keys())

    # first index of letter
    occ = {}

    idx = 0
    for c in letters:
        occ[c] = idx
        idx += A[c]
    del idx, A

    return occ

def calc_checkpoints(s, step):
    """ count the number of letters for each tally matrix step and
        return list of the counts"""
    A = {}  # letter count
    C = []  # checkpoints
    for i, c in enumerate(s):
        if i % step == 0:
            C.append(A.copy())
        if A.get(c):
            A[c] += 1
        else:
            A[c] = 1
    return C


def calc_sa_checkpoints(sa, steps):
    """ count the number of letters for each suffix array step and
            return list of the counts"""
    sa_temp = []
    for i in range(len(sa)):
        if (sa[i] % steps) != 0:
            sa_temp.append(None)
        else:
            sa_temp.append(sa[i])
    return sa_temp


def count_letter_with_checkpoints(C, step, s, idx, letter):
    """ Count the number of a letter upto idx in s using checkpoints.

    Arguments:
    C      -- is the list of checkpoints
    step   -- is the step of the checkpoints
    s      -- the transformed string
    idx    -- count upto this position
    letter -- count for this letter
    """

    # find the nearest checkpoint for idx
    check = int((idx + (step / 2)) / step)
    if check >= len(C):
        check = len(C) - 1
    pos = check * step

    # count of the letter s[idx] upto pos (not included)
    count = C[check].get(letter)
    if count == None:
        count = 0

    # range between pos and idx
    if pos < idx:
        r = range(pos, idx)
    else:
        r = range(idx, pos)

    # count of letters between pos, idx
    k = 0
    for i in r:
        if letter == s[i]:
            k += 1

    # calculate the letter count upto idx (not included)
    if pos < idx:
        count += k
    else:
        count -= k

    return count