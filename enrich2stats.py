'''
enrich2stats.py - Extract and calculate summaries from annotation enrichment testing output
===========================================================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Calculate summary statistics from annotation enrichment testing input.

Usage
-----

.. Example use case

Example::

   python enrich2stats.py

Type::

   python enrich2stats.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
from collections import *
import pandas as pd
import re


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

    parser.add_option("--method", dest="method", type="choice",
                      choices=["merge", "stat", "append"],
                      help="task to perform, either generate "
                      "summary statistics or merge them")

    parser.add_option("--statistic", dest="stat", type="choice",
                      choices=["median", "percent", "mean", "summary"],
                      help="output statistic from the input data")

    parser.add_option("--data-column", dest="data_col", type="int",
                      help="column number containg data over which"
                      " summary stats will be calculated.  Uses "
                      "0-based indexing.")
                      
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    parser.set_defaults(stat="mean",
                        data_col=1)
                        
    infile = argv[-1]

    if options.method == "stat":
        # assume tab separation and header column
        indf = pd.read_table(infile, sep="\t", header=0,
                             index_col=None)

        if options.stat == "median":
            stat_out = indf.iloc[:, options.data_col].median()
        elif options.stat == "percent":
            max_val = indf.iloc[:, options.data_col].max()
            stat_out = (indf.iloc[:, options.data_col]/max_val) * 100.0
        elif options.stat == "mean":
            stat_out = indf.iloc[:, options.data_col].mean()
        elif options.stat == "summary":
            # output a table of summary stats
            median = indf.iloc[:, options.data_col].median()
            min_val = indf.iloc[:, options.data_col].min()
            max_val = indf.iloc[:, options.data_col].max()
            mean_val = indf.iloc[:, options.data_col].mean()
            first_q, last_q = indf.iloc[:, options.data_col].quantile([0.25,
                                                                       0.75])
            stat_out = pd.Series({"Median": median,
                                  "Mean": mean_val,
                                  "Min": min_val,
                                  "Max": max_val,
                                  "25th quantile": first_q,
                                  "75th quantile": last_q})

        if type(stat_out) != float:
            stat_out.to_csv(options.stdout, sep="\t", index_label="Summary")
        else:
            options.stdout.write("{}\n".format(stat_out))

    elif options.method == "merge":
        # parse each file, given the appropriate expectation
        infiles = infile.split(",")
        pval_file = [pf for pf in infiles if re.search("pval", pf)][0]
        with IOTools.openFile(pval_file, "r") as pfile:
            p_val = float(pfile.read())

        boundary_file = [bf for bf in infiles if re.search("boundary",
                                                           bf)][0]
        with IOTools.openFile(boundary_file, "r") as bfile:
            boundary = float(bfile.read())

        summary_file = [sf for sf in infiles if re.search("summary",
                                                          sf)][0]

        # the stats table may be empty if there was no enrichment??
        try:
            stat_df = pd.read_table(summary_file, sep="\t", header=None,
                                    index_col=0).T
        except ValueError:
            # need to convert to a dataframe for output to stdout
            stat_df = pd.Series({"Median": 0,
                                 "Mean": 0,
                                 "Min": 0,
                                 "Max": 0,
                                 "25th quantile": 0,
                                 "75th quantile": 0}).T

        stat_df["P"] = p_val
        stat_df["boundary_size"] = boundary

        # extract the annotation name from the input file and extra
        # columns - useful for multi indexing later

        name = summary_file.split("/")[-1].split(".")[0].split("_")[-1]
        cell_type = name.split("-")[0]
        annotation = name.split("-")[-1]

        stat_df["cell_type"] = cell_type
        stat_df["annotation"] = annotation

        stat_df.to_csv(options.stdout, sep="\t",
                       index=None)

    elif options.method == "append":
        file_list = infile.split(",")
        file_0 = file_list.pop(0)
        E.info("Appending data from file {}".format(file_0))

        df = pd.read_table(file_0,
                           sep="\t",
                           header=0,
                           index_col=None)
        for mfile in file_list:
            E.info("Appending data from file {}".format(mfile))
            _df = pd.read_table(mfile, sep="\t",
                                header=0, index_col=None)
            df = df.append(_df)

        df.to_csv(options.stdout, sep="\t", index=None)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
