'''
snps2snpset.py - Get SNP co-ordinates from a SNP list and .bim file
===================================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Get SNP co-ordinates from a .bim file for a list of SNPs

Usage
-----

.. Example use case

Example::

   python snps2snpset.py

Type::

   python snps2snpset.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
from collections import *
import pandas as pd


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

    parser.add_option("--bim-file", dest="bim_file", type="string",
                      help="Plink .bim file containing SNP positions"
                      "If multiple files are provided (comma separated)"
                      "then all SNPs in all files will be used.")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]
    snp_pos = OrderedDict()

    E.info("Parsing SNP list")

    with IOTools.openFile(infile, "r") as snpfile:
        snpset = snpfile.readlines()

    snpset = [sx.rstrip("\n") for sx in snpset]

    E.info("SNPs in SNP set: %i" % len(snpset))
    E.info("Parsing SNP position information")
    if len(options.bim_file.split(",")) > 1:
        E.info("multiple SNP map files detected")
        for bfile in options.bim_file.split(","):
            with IOTools.openFile(bfile, "r") as ofile:
                for line in ofile:
                    chrom = "chr" + line.split("\t")[0]
                    start = line.split("\t")[3]
                    snp = line.split("\t")[1]
                    try:
                        snp_pos["SNP"].append(snp)
                        snp_pos["Chrom"].append(chrom)
                        snp_pos["BP"].append(start)
                    except KeyError:
                        snp_pos["SNP"] = [snp]
                        snp_pos["Chrom"] = [chrom]
                        snp_pos["BP"] = [start]

    else:
        with IOTools.openFile(options.bim_file, "r") as ofile:
            for line in ofile:
                chrom = "chr" + line.split("\t")[0]
                start = line.split("\t")[3]
                snp = line.split("\t")[1]
                try:
                    snp_pos["SNP"].append(snp)
                    snp_pos["Chrom"].append(chrom)
                    snp_pos["BP"].append(start)
                except KeyError:
                    snp_pos["SNP"] = [snp]
                    snp_pos["Chrom"] = [chrom]
                    snp_pos["BP"] = [start]

    all_df = pd.DataFrame(snp_pos)
    all_df.index = all_df["SNP"]

    snp_df = all_df.loc[snpset]
    snp_df.dropna(axis=0, how='any', inplace=True)
    snp_df.to_csv(options.stdout, sep="\t", index=None)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
