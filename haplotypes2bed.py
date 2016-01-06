'''
haplotypes2bed.py - convert blocks of haplotypes to interval format
===================================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Convert haplotype blocks into interval format (BED4)

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

    parser.add_option("--contig-column", dest="contig_col", type="int",
                      help="file column containing contig information."
                      "Assume 1-based index")

    parser.add_option("--start-column", dest="start_col", type="int",
                      help="file column containing haplotype start "
                      "position. Assume 1-based index")

    parser.add_option("--end-column", dest="end_col", type="int",
                      help="file column containing haplotype end "
                      "position. Assume 1-based index")

    parser.add_option("--snp-column", dest="snp_col", type="int",
                      help="file column containing SNPs within the "
                      "haplotype. Assume 1-based index")

    parser.add_option("--header", dest="header", action="store_true",
                      help="If input file has a header column")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]

    with IOTools.openFile(infile, "r") as ofile:
        lcount = 0
        for line in ofile:
            if options.header and lcount == 0:
                pass
                lcount += 1
            else:
                contig = line.split("\t")[options.contig_col - 1]
                start = line.split("\t")[options.start_col - 1]
                end = line.split("\t")[options.end_col - 1]
                snps = line.split("\t")[options.snp_col - 1].rstrip("\n")
                snps = snps.replace("|", ";")

                options.stdout.write("chr%s\t%s\t%s\t%s\n" % (contig,
                                                              start,
                                                              end,
                                                              snps))
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
