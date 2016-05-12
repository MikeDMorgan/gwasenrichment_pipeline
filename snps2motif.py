'''
snps2motif.py - test SNPs for motif disrupting capability with motifbreakR
===================================================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This tool takes an input list of SNPs or file with a column of SNP ids
and tests for their ability to disrupt known sequence/transcription
factor motifs.  It uses the bioconductor package motifbreakR.

Usage
-----

.. Example use case

Example::

   python snps2motif.py

Type::

   python snps2motif.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri as py2ri
import rpy2.rinterface as ri
import rpy2.robjects as ro
import pandas as pd


def testMotifsDisruption(snp_list, save_path, scripts_dir,
                         r_script, motif_dbname="JASPAR",
                         species="Hsapiens", p_threshold=0.001,
                         motif_pwm=None):
    '''
    Test a list of SNPs for their disruptive influence on known
    DNA binding factor motifs.

    Arguments
    ---------
    snp_list: list
      SNP ids to test for disruptive effects

    save_path: string
      PATH to save image files to

    r_script: string
      R script with motifbreakR workflow

    scripts_dir: string
      PATH to directory containing R analysis scripts

    motif_dbname: string
      A valid motif DB name recognised by motifbreakR

    species: string
      Species to screen motifs from.  If None motifs for
      all species are used.

    p_threshold: float
      p-value threshold for reporting disrupted motifs

    motif_pwm: string
      PATH to file containing positional weight matrices
      of additional motifs to analyse

    Returns
    -------
    motifs_df: pandas.Core.DataFrame
      dataframe of SNPs disrupting specific motifs, with
      degree of disruption motif location
    '''

    # put snp ids into an S4 vector object
    # don't need to convert other arguments as
    # they are standard data type, i.e. int, string, etc

    py2ri.activate()

    # rpy2 doesn't like Python set objects, convert
    # to list first
    r_snps = ro.StrVector([x for x in set(snp_list)])
    R.assign("snp.ids", r_snps)

    E.info("R logging information can be found "
           "in motifbreakR.log")
    # capture output messages in sink file
    R('''sink(file="motifbreakR.log", append=T)''')
    R('''source("{}/{}")'''.format(scripts_dir,
                                   r_script))

    if motif_pwm:
        E.info("Adding motif {} to DB".format(motif_pwm))
        add_motif = ro.vectors.BoolVector([True])
    else:
        add_motif = ro.vectors.BoolVector([False])

    R.assign("add.motif", add_motif)

    E.info("Executing motifbreakR R script: {}".format(r_script))
    R('''resdf <- SnpOnMotif(snp_ids=snp.ids, r_scripts="%(scripts_dir)s",'''
      '''output_dir="%(save_path)s",'''
      '''add_motif=add.motif, '''
      '''motif_pwm="%(motif_pwm)s")''' % locals())
    R('''sink(file=NULL)''')
    res_df = py2ri.ri2py(R["resdf"])

    return res_df


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

    parser.add_option("--snp-column", dest="snp_col", type="int",
                      help="column containing SNP ids.  If this is "
                      "only a list of SNPs, then this is column 1")

    parser.add_option("--header", dest="header", action="store_true",
                      help="Does input file contain a header?")

    parser.add_option("--R-scripts-directory", dest="scripts_dir",
                      type="string", help="directory containing R "
                      "script to run motifbreakR analysis")

    parser.add_option("--R-script", dest="r_script", type="string",
                      help="R script with motifbreakR analysis "
                      "workflow")

    parser.add_option("--additional-motif", dest="add_motif",
                      type="string", help="file containing PWM(s) "
                      "of additional motifs to analyse for disruption")

    parser.add_option("--image-directory", dest="image_dir", type="string",
                      help="directory to save plots in")


    parser.set_defaults(snp_col=1,
                        header=False,
                        add_motif=None)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)


    infile = argv[-1]
    if options.header:
        indf = pd.read_table(infile, sep="\t", index_col=None,
                             header=0)
    else:
        indf = pd.read_table(infile, sep="\t", index_col=None,
                             header=None)

    if type(indf) == pd.core.series:
        snp_ids = indf.values
    else:
        snp_ids = indf.iloc[:, options.snp_col].values

    out_df = testMotifsDisruption(snp_list=snp_ids,
                                  save_path=options.image_dir,
                                  scripts_dir=options.scripts_dir,
                                  r_script=options.r_script,
                                  motif_pwm=options.add_motif)


    out_df.to_csv(options.stdout, sep="\t",
                  index=None)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
