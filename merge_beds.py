'''
merge_beds.py - merge multiple bed files together
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. merge multiple bed files together, use file names as annotation names

Usage
-----

.. Example use case

Example::

   python merge_beds.py

Type::

   python merge_beds.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import re
import CGAT.IOTools as IOTools
import CGAT.Bed as Bed
import itertools


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

    parser.add_option("--regex-filename", dest="re_name", type="string",
                      help="regex for filename component to be used "
                      "as the annotation label")

    parser.set_defaults(
        re_name="(.+)",
    )

    (options, args) = E.Start(parser, argv=argv)

    infiles = argv[-1]
    bed_files = infiles.split(",")

    if len(bed_files) == 1:
        raise IOError("Only one file detected, cannot merge "
        "a single bed file")
    else:
        pass

    # get the regex for annotation names,
    rx = re.compile(options.re_name)
    annot_names = [[y for y in rx.search(x).groups()] for x in bed_files]
    annot_names = list(itertools.chain(*annot_names))

    # output as BED4 format: chr, start, end, name
    for fx in range(len(bed_files)):
        bfile = bed_files[fx]
        with IOTools.openFile(bfile, "r") as ofile:
            intervals = Bed.iterator(ofile)
            track_name = [bx for bx in annot_names if re.search(bx, bfile)][0]
            for entry in intervals:
                entry.__setitem__("name", track_name)
                options.stdout.write("%s\t%s\t%s\t%s\n" % (entry.contig,
                                                           entry.start,
                                                           entry.end,
                                                           entry.name))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
