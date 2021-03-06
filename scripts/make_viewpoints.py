#!/usr/bin/env python

# HiC-Pro
# Copyleft 2015 Institut Curie
# Author(s): Nicolas Servant, Eric Viara
# Contact: nicolas.servant@curie.fr
# This software is distributed without any guarantee under the terms of the
# GNU General
# Public License, either Version 2, June 1991 or Version 3, June 2007.

"""
Script to keep only valid 3C products - DE and SC are removed
"""

import getopt
import sys
import os
import re
from bx.intervals.intersection import Intersecter, Interval


def usage():
    """Usage function"""
    print("Usage : python mapped_2hic_fragments.py")
    print("-i/--validPairs <list of valid Hi-C pairs as generated by HiC-Pro>")
    print("-f/--fragmentFile <Restriction fragment file (BED)>")
    print("-t/--target <Target file for viewpoints (BED)>")
    print("[-e/--exclusion] <Size of the regions upstream/downstream the capture sites to discard>")
    print("[-o/--output] <Output file. Default is stdout>")
    print("[-v/--verbose] <Verbose>")
    print("[-h/--help] <Help>")
    return

class Error_handling:


def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "i:f:t:e:o:vh",
            ["validPairsFile=",
             "fragmentsFile=",
             "targetFile=",
             "exclusionSize=",
             "output=", 
             "verbose", "help"])
    except getopt.GetoptError, err:
        print("GetoptError: {}\n".format(str(err))
        usage()
        sys.exit(-1)
    return opts


def load_BED(in_file, exclusionSize=0, verbose=False):
    """
    Read a BED file and store the intervals in a tree

    Intervals are zero-based objects. The output object is a hash table with
    one search tree per chromosome

    in_file = input file [character]
    verbose = verbose mode [logical]

    """
    x = {}
    x_ex = {}
    if verbose:
        print("## Loading BED file {} ...".format(in_file), file=sys.stderr)
    nline = 0
    with open(in_file) as bed_handle:
        for line in bed_handle:
            if nline % 1000000 == 0 and verbose:
                print("{} million lines loaded".format(int(nline/1000000))
            nline += 1
            bedtab = line.split("\t")
            try:
                chromosome, start, end, name = bedtab[:4]
            except ValueError:
                print("Warning : wrong input format in line {}. Not a BED file !?".format(nline), file=sys.stderr)
                continue
            
            # BED files are zero-based as Intervals objects
            start = int(start)  # + 1
            end = int(end)
            name = name.strip()
            if chromosome in x:
                tree = x[chromosome]
                tree.add_interval(Interval(start, end, value={'name': name}))
            else:
                tree = Intersecter()
                tree.add_interval(Interval(start, end, value={'name': name}))
                x[chromosome] = tree             
            ## Exclusion regions
            if exclusionSize > 0:
                if chromosome in x_ex:
                    tree_ex = x_ex[chromosome]
                    tree_ex.add_interval(Interval(start - int(exclusionSize), start, value={'name': str(name) + "_up"}))
                    tree_ex.add_interval(Interval(end, end + int(exclusionSize), value={'name': str(name) + "_dwn"}))
                else:
                    tree_ex = Intersecter()
                    tree_ex.add_interval(Interval(start - int(exclusionSize), start, value={'name': str(name) + "_up"}))
                    tree_ex.add_interval(Interval(end, end + int(exclusionSize), value={'name': str(name) + "_dwn"}))
                    x_ex[chromosome] = tree_ex             
    return (x, x_ex)
    
    

def get_overlapping_fragment(frag, chrom, pos, quiet=False):
    """
    Intersect a given read with the set of restriction fragments

    ##
    frag = the fragments [hash]
    chrom = the chromosome to look at [character]
    read = the read to intersect [AlignedRead]

    """
    if chrom in frag:
        # Overlap with the start of the read (zero-based)
        ifrag = frag[chrom].find(int(pos), int(pos+1))
        if len(ifrag) > 1:
            if not quiet:
                print(file=sys.stderr, "Warning : {} fragments found for read at {}:\
                                        {} -skipped {}".format(len(ifrag), chrom, pos, ifrag))
            return None
        elif len(ifrag) == 0:
            if not quiet:
                print(file=sys.stderr, "Warning - no fragments found for read at {}:\
                                        {} -skipped".format(chrom, pos))
            return None
        else:
            return ifrag[0]
    else:
        if not quiet:
            print(file=sys.stderr, "Warning - no fragments found for read at {}:\
                                    {} -skipped".format(chrom, pos))
        return None


if __name__ == "__main__":
    # Read command line arguments
    opts = get_args()
    verbose = False
    output = None
    exclusionSize = 0

    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--validPairs"):
            validPairsFile = arg
        elif opt in ("-f", "--fragmentFile"):
            fragmentFile = arg
        elif opt in ("-t", "--targetFile"):
            targetFile = arg
        elif opt in ("-e", "--exclusion"):
            exclusionSize = arg
        elif opt in ("-o", "--output"):
            output = arg
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # Verbose mode
    if verbose:
        print(file=sys.stderr, "## make_viewpoints.py")
        print(file=sys.stderr, "## validPairsFile={}".format(validPairsFile))
        print(file=sys.stderr, "## fragmentFile={}".format(fragmentFile))
        print(file=sys.stderr, "## targetFile={}".format(targetFile))
        print(file=sys.stderr, "## exclusionSize={}".format(exclusionSize))
        print(file=sys.stderr, "## verbose={}\n".format(verbose))


    # Read the BED files
    if verbose:
        print(file=sys.stderr, "## Loading data ...")

    resFrag = load_BED(fragmentFile, verbose)[0]
    (target, exclu) =  load_BED(targetFile, exclusionSize=exclusionSize, verbose=verbose)

    # Read the validPairs file
    if verbose:
        print(file=sys.stderr, "## Opening file {} ...".format(validPairsFile))
   
    nline = 0
    repdict = {}

    c_c_counter, r_r_counter, c_r_counter, ua_counter, exclu_counter = 0, 0, 0, 0, 0

    with open(validPairsFile) as in_handle:
        if nline % 1000000 == 0 and verbose:
            print(file=sys.stderr, "{} million lines processed".format(int(nline/1000000)))
        for line in in_handle:
            nline += 1
            intab = line.split("\t")
            try:
                readname, r1_chr, r1_start, r1_strand, r2_chr, r2_start, r2_strand = intab[:7]
            except ValueError:
                print(file=sys.stderr, "Warning : wrong input format in line {}\
                                        .Not a validPairs file !?".format(nline))
                continue

            r1_resfrag, r2_resfrag, reporter = None, None, None
            
            ## Intersect with target
            if len(exclu) > 0:
                r1_resfrag = get_overlapping_fragment(exclu, r1_chr, int(r1_start), quiet=True)
                r2_resfrag = get_overlapping_fragment(exclu, r2_chr, int(r2_start), quiet=True)
                
            if r1_resfrag is None and r2_resfrag is None:
                r1_resfrag = get_overlapping_fragment(target, r1_chr, int(r1_start), quiet=True)
                r2_resfrag = get_overlapping_fragment(target, r2_chr, int(r2_start), quiet=True)

                if r1_resfrag is not None:
                    capture = r1_resfrag  
                else:
                    r1_resfrag = get_overlapping_fragment(resFrag, r1_chr, int(r1_start))
                    if  r1_resfrag is not None:
                        reporter = r1_resfrag
                    else:
                        ua_counter += 1
                        #print >> sys.stderr, "Warning : reads [", r1_chr, ":", r1_start, "] do not overlap with a capture nor a restriction fragment !"

                
                if r2_resfrag is not None:
                    capture = r2_resfrag
                else:
                    r2_resfrag = get_overlapping_fragment(resFrag, r2_chr, int(r2_start))
                    if r2_resfrag is not None:
                        reporter = r2_resfrag
                    else:
                        ua_counter += 1
                        #print >> sys.stderr, "Warning : reads [", r2_chr, ":", r2_start, "] do not overlap with a capture nor a restriction fragment !"
            else:
                exclu_counter += 1
        
            ## Counts
            if capture is not None and reporter is not None: 
                c_r_counter += 1
                if not capture.value['name'] in repdict:
                    repdict[capture.value['name']] = {}                                                                                                                     

                if reporter.value['name'] in repdict[capture.value['name']]:
                    repdict[capture.value['name']][reporter.value['name']]['count'] += 1                                                                             
                else:
                    repdict[capture.value['name']][reporter.value['name']] = {'chr':r1_chr, 'start':reporter.start, 'end':reporter.end, 'count':1}                  
            elif capture is not None and reporter is None:
                c_c_counter += 1
            elif capture is None and reporter is not None:
                r_r_counter += 1

    in_handle.close()

    ## Write
    if output is not None:
        sys.stdout = open(output, 'w')

    for k in repdict:
        print("track type=bedGraph name='hicpro {k}' description='hicpro {k}' visibility=full\
               color=200,100,0 altColor=0,100,200 priority=20".format(k, k))
        for key, value in repdict[k].iteritems():
            print("{}\t{}\t{}\t{}".format(value['chr'], str(value['start']),\
                                          str(value['end']), str(value['count'])))


    ## stats
    print(file=sys.stderr, "CAP-REP read pairs = {}".format(c_r_counter))
    print(file=sys.stderr, "CAP-CAP read pairs = {}".format(c_c_counter))
    print(file=sys.stderr, "REP-REP read pairs = ".format(r_r_counter))
    print(file=sys.stderr, "Excluded reads =".format(exclu_counter))
    print(file=sys.stderr, "UA reads =".format(ua_counter))
