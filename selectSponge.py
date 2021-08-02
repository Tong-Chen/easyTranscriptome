#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from __future__ import division, with_statement
'''
Copyright 2015, 陈同 (chentong_biology@163.com).  
===========================================================
'''
__author__ = 'chentong & ct586[9]'
__author_email__ = 'chentong_biology@163.com'
#=========================================================
desc = '''
Program description:

Sponge file:

    miR-18a-3p      ENST00000246190 1       Coding
    miR-18a-3p      ehbio_trans.530.1       1       lncRNA
    miR-18a-3p      ehbio_trans.530.3       1       lncRNA
    miR-18a-3p      ehbio_trans.530.5       1       lncRNA
    miR-18a-3p      ehbio_trans.300.1       1       lncRNA
    miR-18a-3p      ehbio_trans.573.1       1       lncRNA
    miR-18a-3p      ehbio_trans.22.4        1       lncRNA
    miR-18a-3p      ENST00000484569 1       lncRNA
    miR-18a-3p      ENST00000465000 1       lncRNA
    miR-18a-5p      ENST00000338244 1       Coding

Correlation file:



'''

import sys
import os
from json import dumps as json_dumps
from time import localtime, strftime 
timeformat = "%Y-%m-%d %H:%M:%S"
from optparse import OptionParser as OP
#from multiprocessing.dummy import Pool as ThreadPool

#from bs4 import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')


debug = 0


def cmdparameter(argv):
    if len(argv) == 1:
        global desc
        print(desc)
        cmd = 'python ' + argv[0] + ' -h'
        os.system(cmd)
        sys.exit(1)
    usages = "%prog -i file"
    parser = OP(usage=usages)
    parser.add_option("-i", "--sponge-file", dest="filein",
        metavar="FILEIN", help="")
    parser.add_option("-c", "--correlation-file", dest="correlation",
        metavar="FILEIN", help="")
    parser.add_option("-o", "--output-prefix", dest="op",
        metavar="FILEIN", help="")
    parser.add_option("-v", "--verbose", dest="verbose",
        action="store_true", help="Show process information")
    parser.add_option("-d", "--debug", dest="debug",
        default=False, action="store_true", help="Debug the program")
    (options, args) = parser.parse_args(argv[1:])
    assert options.filein != None, "A filename needed for -i"
    return (options, args)
#--------------------------------------------------------------------

def readInSponge(sponge):
    spongeD = {}
    for line in open(sponge):
        miR, gene, cnt, type1 = line.strip().split()
        if miR not in spongeD:
            spongeD[miR] = {}
        if type1 not in spongeD[miR]:
            spongeD[miR][type1] = []
        spongeD[miR][type1].append(gene)
    return spongeD
#----------------------------------------------
def readInCor(correlation):
    correlationD= {}
    header = 1
    for line in open(correlation):
        lineL = line.strip().split()
        gene1,  gene2 = lineL[:2]
        correlationD[(gene1, gene2)] = lineL[2]
        correlationD[(gene2, gene1)] = lineL[2]
    return correlationD
#--------------------------------------------

def main():
    options, args = cmdparameter(sys.argv)
    #-----------------------------------
    sponge = options.filein
    correlation = options.correlation
    op = options.op
    verbose = options.verbose
    global debug
    debug = options.debug
    #-----------------------------------
    spongeD = readInSponge(sponge)
    correlationD = readInCor(correlation)
    
    network = op + '.network.txt'
    network_fh = open(network, 'w')
    print("\t".join(['Source', "Target", 'Correlation']), file=network_fh)
    table   = op + '.table.txt'
    table_fh = open(table, 'w')
    print("Gene\tClass", file=table_fh)

    for miR, typeD in spongeD.items():
        keyL = list(typeD.keys())
        assert len(keyL)==2, "Unsupported data" 
        key1 = keyL[0]
        key2 = keyL[1]
        geneL1 = typeD[key1]
        geneL2 = typeD[key2]
        for gene1 in geneL1:
            for gene2 in geneL2:
                if (gene1, gene2) in correlationD:
                    print('\t'.join([miR, gene1, '0']), file=network_fh)
                    print('\t'.join([miR, gene2,'0']), file=network_fh)
                    print('\t'.join([gene1, gene2, correlationD[(gene1, gene2)]]), file=network_fh)
                    print('\t'.join([gene1,key1]), file=table_fh)
                    print('\t'.join([gene2,key2]), file=table_fh)


    network_fh.close()
    table_fh.close()
    #-----------end close fh-----------
    ###--------multi-process------------------
    #pool = ThreadPool(5) # 5 represents thread_num
    #result = pool.map(func, iterable_object)
    #pool.close()
    #pool.join()
    ###--------multi-process------------------

if __name__ == '__main__':
    startTime = strftime(timeformat, localtime())
    main()
    endTime = strftime(timeformat, localtime())
    ###---------profile the program---------
    #import profile
    #profile_output = sys.argv[0]+".prof.txt")
    #profile.run("main()", profile_output)
    #import pstats
    #p = pstats.Stats(profile_output)
    #p.sort_stats("time").print_stats()
    ###---------profile the program---------

