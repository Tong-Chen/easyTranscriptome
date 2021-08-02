#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, with_statement
'''
Copyright 2013, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Functional description:
    This is designed to extract FASTA seq by names.
'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

debug = 0

def fprint(content):
    print(json_dumps(content,indent=1))

def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print(desc, file=sys.stderr)
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--input-file", dest="filein",
        metavar="FILEIN", help="The FASTA file")
    parser.add_option("-s", "--separator", dest="sep",
        default="IdoNotThinkThisWillAppear, DoyouThinkSo", 
        metavar="SEPARATOR", help="The separator used to get ID names. \
Default full line (no splitting) except leading > and trailing '\\n' is \
used as ID names. Please use <tab> to specify '\\t' as separtor.")
    parser.add_option("-F", "--first-x-words", dest="count",
        default=1, help="Default 1 means extracting the first \
word before separator. Accept other number (x) to extract the \
first x words.")
    parser.add_option("-n", "--name-list", dest="name",
        help="One or several columns file containing ID lists in one column.")
    parser.add_option("-c", "--name-col-index", dest="name_col_ix",
        default=1, type='int', 
        help="Specify the columns containing IDs. Default 1 representing the first column.")
    parser.add_option("-r", "--rename-col-index", dest="rename_col_ix",
        default=0, type='int', 
        help="Specify the columns containing IDs. Default 0 representing no rename.")
    parser.add_option("-v", "--verbose", dest="verbose",
        default=0, help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------


def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    file = options.filein
    global debug
    debug = options.debug
    sep = options.sep
    #print "*%s*" % sep
    if sep == 'tab':
        sep = "\t"
    count = int(options.count)
    nameF = options.name
    nameC = options.name_col_ix-1
    rename_col_ix = options.rename_col_ix
    if rename_col_ix:
        rename_col_ix = rename_col_ix - 1
    else:
        rename_col_ix = nameC
    nameD = dict([[line.strip().split()[nameC], line.strip().split()[rename_col_ix]] \
            for line in open(nameF)])

    if debug:
        print(list(nameD.keys()), file=sys.stderr)
    #print nameD.keys()[:5]
    verbose = options.verbose
    #-----------------------------------
    if file == '-':
        fh = sys.stdin
    else:
        fh = open(file)
    #--------------------------------
    output = 0
    tmpL = []
    oldline = ''
    for line in fh:
        #print(line)
        if line[0] == '>':
            if output and tmpL and oldline:
                print(oldline.strip())
                print(''.join(tmpL))
                tmpL = []
                output = 0
            key = sep.join(line[1:].strip().split(sep)[:count])
            #print key
            #break
            oldline = line
            #if key in nameD:
            #   nameD.pop(key)
            #    output = 1
            output = nameD.pop(key, 0)
            if output and rename_col_ix != nameC:
                oldline = oldline.strip() + ' ' + output
            if debug:
                print(key, output, file=sys.stderr)
        else:
            if output:
                tmpL.append(line.strip())
    #-------------END reading file----------
    if output and tmpL and oldline:
        print(oldline, end=' ')
        print(''.join(tmpL))
        tmpL = []
        output = 0
    #----close file handle for files-----
    if file != '-':
        fh.close()
    #-----------end close fh-----------
    if nameD:
        print("The following IDs have no sequences found", file=sys.stderr)
        print('\n'.join(list(nameD.keys())), file=sys.stderr)
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------
    if verbose:
        print("--Successful %s" % strftime(timeformat, localtime()), file=sys.stderr)

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    fh = open('python.log', 'a')
    print("%s\n\tRun time : %s - %s " % \
        (' '.join(sys.argv), startTime, endTime), file=fh)
    fh.close()
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------


