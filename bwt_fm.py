import collections
import time, psutil, os
from operator import itemgetter
###########################################################################################
# Global variables
pattern1 = 'TATATAA'
pattern2 = 'CCGGCTAT'
pattern3 = 'TTCACTACTCTCA'

# Canis lupus familiaris genome, chromosome 1
file_name = '4615_ref_ASM154086v1_chr25.fa'

###########################################################################################
argsort = lambda l: [i for i, _ in sorted(enumerate(l), key=itemgetter(1))]

def diff_time(start, end):
    return float((end - start) * 1000)

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

class BurrowsWheeler():
    EOS = "\0"

    # EOS = "#" # a visible end marker
    def suffix_array(self, arr):
        arr_size = len(arr)
        arr_int = {v: k for k, v in enumerate(sorted(set(arr)))}
        arr = [arr_int[x] for x in arr]
        arr.append(-1)
        suf = [[i, arr[i], arr[i + 1]] for i in range(arr_size)]
        # suf = radixSort(suf)
        suf.sort(key=itemgetter(1, 2))
        idx = [0] * arr_size
        k = 2
        while k < arr_size:
            r = 0
            prev_r = suf[0][1]
            for i in range(arr_size):
                if suf[i][1] != prev_r or suf[i - 1][2] != suf[i][2]:
                    r += 1
                prev_r = suf[i][1]
                suf[i][1] = r
                idx[suf[i][0]] = i
            for i in range(arr_size):
                next_idx = suf[i][0] + k
                suf[i][2] = suf[idx[next_idx]][1] if next_idx < arr_size else -1
            suf.sort(key=itemgetter(1, 2))
           # suf = radixSort(suf)
            k <<= 1
        return [x[0] for x in suf]

    def bwt(self,input_string):
        data_ref = self.suffix_array(input_string)
        #self.sa = data_ref
        return (x - 1 for x in data_ref), data_ref.index(0), data_ref


def calc_checkpoints(s, step):
    """ count the number of letters for each step and
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

class SuffixArrayBurrowsWheeler(BurrowsWheeler):

    def transform(self, s):
        """ Burrow-Wheeler transform with SuffixArray,
            similar to SuffixTree implementations. """
        assert self.EOS not in s, "Input string cannot contain null character (%s)" % self.EOS

        # add end of text marker
        s += self.EOS

        # table of suffixes
        rotations = [s[i:] for i in range(len(s))]

        sa = {}
        for i, q in enumerate(rotations):
            sa[q] = i
        od = collections.OrderedDict(sorted(sa.items()))

        # sort the suffixes
        rotations.sort()

        # get the length of ordered suffixes
        k = len(rotations)

        r = [0] * k
        for i in range(k):
            l = len(rotations[i])
            if l == k:
                r[i] = self.EOS
            else:
                r[i] = s[-l - 1]
        r = ''.join(r)

        return r, list(od.values())

class FMSimpleIndex(object):
    def __init__(self, data):
        self.data = bw.transform(data)
        self.offset = {}
        self._build(data)

    def _build(self, data):
        """ build the index """
        self.occ = calc_first_occ(self.data)

    def _occ(self, qc):
        """ get the first occurance of letter qc in left-column"""
        c = self.occ.get(qc)
        if c == None:
            return 0
        return c

    def _count(self, idx, qc):
        """ count the occurances of letter qc (rank of qc) upto position idx """
        if not qc in self.occ.keys(): return 0
        c = 0
        for i in range(idx):
            if self.data[i] == qc:
                c += 1
        return c

    def _lf(self, idx, qc):
        """ get the nearset lf mapping for letter qc at position idx """
        o = self._occ(qc)
        c = self._count(idx, qc)
        return o + c

    def _walk(self, idx):
        """ find the offset in position idx of transformed string
            from the beginning """

        # walk to the beginning using lf mapping
        # this is same as inverse of burrow wheeler transformation
        # from arbitrary location
        r = 0
        i = idx
        while self.data[i] != bw.EOS:
            if self.offset.get(i):
                # we have cached the location and can use it
                r += self.offset[i]
                break
            r += 1
            i = self._lf(i, self.data[i])

        # save the offset of some idx for faster searches
        if not self.offset.get(idx):
            self.offset[i] = r
        return r

    def suffix(self, i):
        count = 0
        while self.sa[i] == None:
            i = self._lf(i, self.data[i])
            count = count + 1
        return self.sa[i] + count

    def bounds(self, q):
        """ find the first and last suffix positions for query q """
        top = 0
        bot = len(self.data)
        for i, qc in enumerate(q[::-1]):
            top = self._lf(top, qc)
            bot = self._lf(bot, qc)
            if top == bot: return (-1, -1)
        return (top, bot)

    def search(self, q):
        """ search the positions of query q """

        # find the suffixes for the query
        top, bot = self.bounds(q)
        matches = []
        # find the location of the suffixes
        # by walking the reverse text from that position
        # with lf mapping
        for i in range(top, bot):
            # pos = self._walk(i)
            pos = self.suffix(i)
            matches.append(pos)
        return sorted(matches)

    def count(self, q):
        """ count occurances of q in the index """
        top, bot = self.bounds(q)
        return bot - top

class FMCheckpointing(FMSimpleIndex):
    """ creates LF index with checkpoints """

    def __init__(self, data):
        bwt_ref, idx, self.sa = bw.bwt(data)
        encoded = "".join(data[x] for x in bwt_ref)
        self.data = encoded
        self.offset = {}
        self.step = tally_step
        self.step_sa = sa_step
        self._build()

    def _build(self):
        """ build the index """
        self.occ = calc_first_occ(self.data)
        self.C = calc_checkpoints(self.data, self.step)
        self.sa = calc_sa_checkpoints(self.sa, self.step_sa)

    def _count(self, idx, qc):
        """ count the occurances of letter qc (rank of qc) upto position idx """
        count = count_letter_with_checkpoints(self.C, self.step, self.data, idx, qc)
        return count

def prepare_file(file_name):
    global data

    print('\nPreparing file ' + file_name + '...')
    file = open(file_name, 'r')
    file.readline()
    data = file.read().replace('\n', '')
    data = data + '#'
    file.close()
    print('Prepared!')

def file_processing(file_name):
    global bw

    prepare_file(file_name)

    pid = os.getpid()
    ps = psutil.Process(pid)

    t_start_bwt = time.clock()
    bw = SuffixArrayBurrowsWheeler()

    idx = FMCheckpointing(data)
    t_end_bwt = time.clock()
    t_start_pattern1 = time.clock()
    m1 = idx.search(pattern1)
    t_end_pattern1 = time.clock()

    t_start_pattern2 = time.clock()
    m2 = idx.search(pattern2)
    t_end_pattern2 = time.clock()

    t_start_pattern3 = time.clock()
    m3 = idx.search(pattern3)
    t_end_pattern3 = time.clock()
    print(m1)
    print(m2)
    print(m3)

    print("BWT: %sms" % diff_time(t_start_bwt, t_end_bwt))
    print("FM for pattern %s for tally_step %s and sa_step %s: %sms" % (pattern1, tally_step, sa_step, diff_time(t_start_pattern1, t_end_pattern1)))
    print("FM for pattern %s for tally_step %s and sa_step %s: %sms" % (pattern2, tally_step, sa_step, diff_time(t_start_pattern2, t_end_pattern2)))
    print("FM for pattern %s for tally_step %s and sa_step %s: %sms" % (pattern3, tally_step, sa_step, diff_time(t_start_pattern3, t_end_pattern3)))
    print("Memory usage: {}".format(ps.memory_info()[0] / 2. ** 20))

# Choose string matching algorithm
def main():
    choose_solution()

    file_processing(file_name)

def choose_solution():
    global choice
    global tally_step
    global sa_step
    global structure_name

    print('Choose solution:\n')
    print('1. Non optimized solution\n')
    print('2. Optimized solution\n')

    while True:
        try:
            choice = int(input())
            if choice < 1 or choice > 2:
                raise ValueError

            if choice == 2:
                structure_name = 'optimized solution'
                print('Choose factor for tally matrix: 8, 32, 128 or 512\n')
                tally_step = int(input())
                if tally_step != 8 and tally_step != 32 and tally_step != 128 and tally_step != 512:
                    raise ValueError
                print('Choose factor for suffix array: 4, 16, 64, 254\n')
                sa_step = int(input())
                if sa_step != 4 and sa_step != 16 and sa_step != 64 and sa_step != 254:
                    raise ValueError
            else:
                sa_step = 1
                tally_step = 1

            return choice
        except ValueError:
            print('***ERROR***\n', 'Enter a valid input\n')

main()