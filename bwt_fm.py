import time, psutil, os, sys
import bwt_build
import fm_build
from os.path import isfile
###########################################################################################

""" Default file for processing """
file_name = '4615_ref_ASM154086v1_chr25.fa'

###########################################################################################

""" Calculate difference between start and end time in miliseconds """
def diff_time(start, end):
    return float((end - start) * 1000)

""" Function for reading file passed as argument, if there is one, or default file if no argument is passed """
def prepare_file(file_name):
    if len(sys.argv) == 2:
        if not isfile(sys.argv[1]):
            print("Input file doesn't exist")
            os.abort()
        file_name = sys.argv[1]

    print('\nPreparing file ' + file_name + '...')
    file = open(file_name, 'r')
    file.readline()
    data = file.read().replace('\n', '')
    data = data + '#'
    file.close()
    print('Prepared!')
    return data

""" Core function """
def file_processing(file_name):
    #read file
    data = prepare_file(file_name)

    #read patterns and factors for suffix array and tally matrix
    pattern1, pattern2, pattern3, sa_step, tally_step = input_parameters()

    #memory usage calculation
    pid = os.getpid()
    ps = psutil.Process(pid)

    #start time for bwt creation
    t_start_bwt = time.clock()

    #Create burrows-wheeler transformation, suffix array and tally matrix
    bw = bwt_build.SuffixArrayBurrowsWheeler()
    idx = fm_build.FMCheckpointing(data, bw, sa_step, tally_step)

    #end time for bwt creation
    t_end_bwt = time.clock()

    #calculate time and search for first pattern
    t_start_pattern1 = time.clock()
    m1 = idx.search(pattern1)
    t_end_pattern1 = time.clock()

    # calculate time and search for second pattern
    t_start_pattern2 = time.clock()
    m2 = idx.search(pattern2)
    t_end_pattern2 = time.clock()

    # calculate time and search for third pattern
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

""" Choose string matching algorithm """
def main():
    file_processing(file_name)

""" Function for reading patterns and sa and tally factors """
def input_parameters():

    print('Insert first pattern:\n')
    p1 = input()
    print('Insert second pattern:\n')
    p2 = input()
    print('Insert third pattern:\n')
    p3 = input()
    print('Insert suffix array factor:\n')
    s = input()
    print('Insert tally matrix factor:\n')
    t = input()

    if not p1:
        p1 = 'CTAAATT'
    if not p2:
        p2 = 'GTGTGCCT'
    if not p3:
        p3 = 'AAACCCTAATATG'
    if not s:
        s = 1
    if not t:
        t = 1

    return p1, p2, p3, int(s), int(t)

main()
