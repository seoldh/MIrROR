#!/usr/bin/python

from  lib.LibMIrROR import *
import sys, re, os, datetime, platform, glob, math
from optparse import OptionParser, OptionGroup

###################
# main
###################
if __name__ == "__main__":
    desc = "MIrROR: a tool for metataxonomics with 16S-23S rRNA operon region.\
    INPUTFILE: FASTA/FASTQ/PAF | SAMPLELIST"

    parser = OptionParser(usage='usage: %prog [options]  (-d DBDIR) INPUTFILE', description = desc, add_help_option=False)
    parser.add_option("-d", "--db_dir", dest="databasefolder", action="store",  help="path to 16S-23S rRNA operon database [required]", metavar="DBDIR")
    parser.add_option("-o", "--output_dir", dest="outputfolder", action="store", default='./Result', help="write report to FOLDER (default : ./Result)", metavar="OUTDIR")
    parser.add_option("-t", "--threads", dest="threads", action="store", default=4, type='int',  help="number of threads (default : 4)", metavar="INT")

    mapping_opts = OptionGroup(parser, 'Mapping option', 'this option control for Mapping')
    mapping_opts.add_option("-M", "--minibatch", dest="minibatch", action="store", default="500M", help="the number of query bases loaded to the memory at once. K/M/G/k/m/g suffix is accepted. (default : 500M)", metavar="NUM")
    parser.add_option_group(mapping_opts)

    outformat_opts = OptionGroup(parser, 'Threshold options', 'these options control for Threshold')
    outformat_opts.add_option("-m", "--residuemaches", dest="residuemaches", action="store", default=2500, type='int',  help="number of residue matches in mapping (default : 2500)", metavar="INT")
    outformat_opts.add_option("-b", "--blocklength", dest="blocklength", action="store", default=3500, type='int',  help="alignment block length in mapping (default : 3500)", metavar="INT")
    outformat_opts.add_option("-n", "--Normalization", dest="normalization", action="store_true", help="normalize by 16S-23S rRNA operon copy number")
    parser.add_option_group(outformat_opts)

    graph_opts = OptionGroup(parser, 'Visualization options', 'these options control for Visualization')
    graph_opts.add_option("-K", "--Krona", dest="krona", action="store_true", help="create a Krona plot. Requires KronaTools to be installed")
    graph_opts.add_option("-S", "--Stackedplot", dest="stackedplot", action="store_true", help="create a stacked plot. Requires Python packages to be installed")
    graph_opts.add_option("-V", "--visualized", dest="visualized", action="store_true", help="perform all visualization tasks (Same as: -K -S)")
    parser.add_option_group(graph_opts)

    help_opts = OptionGroup(parser, 'Help')
    help_opts.add_option("-h", "--help", dest="help", action="store_true", help="Show this help message and exit")
    help_opts.add_option("-v", "--version", dest="version", action="store_true", help="Show program's version number and exit\n")
    parser.add_option_group(help_opts)

    (options, args) = parser.parse_args()
    
    if( options.help ):
        parser.print_help()
        sys.exit()
    if( options.version ):
        WriteLog("Info", "%s" % MYVERSION)
        sys.exit()

    if( len(args) == 0 ):
        parser.print_help()
        WriteLog("Error","INPUTFILE required.")
        sys.exit()
    if( options.databasefolder is None):
        WriteLog("Error","Need path to 16S-23S rRNA operon database (-d DBDIR)")
        parser.print_help()
        sys.exit()

    if( not options.databasefolder.startswith("/")):
        options.databasefolder = GETDIR(os.getcwd(),options.databasefolder)                
    if( not os.path.isdir(options.databasefolder) ):
        WriteLog("Error","Not Exist '16S-23S rRNA operon database' folder : %s" % (options.databasefolder))
        parser.print_help()
        sys.exit()

    if(not options.outputfolder.startswith('/') ):
        options.outputfolder = GETDIR(os.getcwd(),options.outputfolder)            
    if( os.path.isdir(os.path.dirname(options.outputfolder)) ):
        if( os.path.isdir(options.outputfolder) ):
            WriteLog("Error","Output folder already exists : %s" % (options.outputfolder))
            parser.print_help()
            sys.exit()
        else:
            os.mkdir(options.outputfolder)
    else:
        WriteLog("Error","Not Found Directory : %s" % (os.path.dirname(options.outputfolder)) )
        parser.print_help()
        sys.exit()

    Main(options,args,sys.argv)
