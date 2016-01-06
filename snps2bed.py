'''
snps2bed.py - merge SNP sets and transform into intervals based on HapMap LD blocks
===================================================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. To merge SNP sets and select representative linkage blocks as
 intervals

Usage
-----

.. Example use case

Example::

   python snps2bed.py

Type::

   python snps2bed.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

def parseDistiLDBlocks(ld_blocks):
    '''
    Parse the DistiLD LD blocks file into intervals
    '''

    interval_dict = {}
    with IOTools.openFile(ld_blocks, "r") as ofile:
        for line in ofile:
            interval = line.split("\t")[0]
            snps = line.split("\t")[1]
            interval_dict[snps] = interval

    return interval_dict


def parsePlinkLDBlocks(ld_blocks):
    '''
    Parse the DistiLD LD blocks file into intervals
    '''

    interval_dict = {}
    with IOTools.openFile(ld_blocks, "r") as ofile:
        for line in ofile:
            contig = line.split("\t")[0]
            start = line.split("\t")[1]
            end = line.split("\t")[2]
            snps = line.split("\t")[-1].rstrip("\n")
            interval = "%s:%s-%s" % (contig, start, end)
            interval_dict[snps] = interval

    return interval_dict


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--ld-blocks", dest="ld_blocks", type="string",
                      help="file containing interval information and "
                      "SNP rs#IDs assigned to each interval")

    parser.add_option("--ld-source", dest="ld", type="choice",
                      choices=["DitiLD", "plink"],
                      help="prorgam used to generate linkage/"
                      "haplotype blocks")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    if options.ld == "DistiLD":
        ld_intervals = parseDistiLDBlocks(options.ld_blocks)
    elif options.ld == "plink":
        ld_intervals = parsePlinkLDBlocks(options.ld_blocks)

    with IOTools.openFile(infile, "r") as snpfile:
        snpset = set([sx.rstrip("\n") for sx in snpfile.readlines()])

    snp_intervals = {}
    # intersect with intervals
    for key in ld_intervals.keys():
        intersect = snpset.intersection(key.split(";"))
        if len(intersect):
            ld_block = ld_intervals[key]
            snp_intervals[ld_block] = intersect
        else:
            pass

    for block in snp_intervals.keys():
        contig = block.split(":")[0]
        start = block.split(":")[-1].split("-")[0]
        end = block.split(":")[-1].split("-")[-1]
        name = ";".join(snp_intervals[block])
        options.stdout.write("%s\t%s\t%s\t%s\n" % (contig,
                                                   start,
                                                   end,
                                                   name))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
