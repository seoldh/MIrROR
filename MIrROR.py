#!/usr/bin/env python3

from  lib.LibMIrROR import *
import sys, os, argparse

###################
# Main
###################
if __name__ == "__main__":
    BANNER = r""" __  __ ___      ____   ___  ____          _   _ 
|  \/  |_ _|_ __|  _ \ / _ \|  _ \  __   _/ | / |
| |\/| || || '__| |_) | | | | |_) | \ \ / / | | |
| |  | || || |  |  _ <| |_| |  _ <   \ V /| |_| |
|_|  |_|___|_|  |_| \_\\___/|_| \_\   \_/ |_(_)_|"""
    
    class MIrRORHelpFormatter(argparse.RawDescriptionHelpFormatter):
        def format_help(self):
            help_text = super().format_help()
            return BANNER + "\n\n" + help_text

    parser = argparse.ArgumentParser(
            usage='%(prog)s [options] (-d DBDIR) INPUTFILE',
            add_help=False,
            epilog="Citation:\nSeol, D. et al. (2022). Microbial Identification Using rRNA Operon Region: Database and tool for metataxonomics with long-read sequence. Microbiol. Spectr., 10(2):e02017-21\nLee,  J. et al. (2026). MIrROR release 02: Expanded and refined 16S-ITS-23S rRNA operon dataset. Sci. Data, 13(1):714",
            formatter_class=MIrRORHelpFormatter
            #formatter_class=argparse.RawDescriptionHelpFormatter
            )
    
    # INPUT
    input_group = parser.add_argument_group('Input')
    input_group.add_argument("inputfiles", nargs='*', metavar="INPUTFILE",
                             help="FASTA/FASTQ/PAF file(s) or a sample list [required]")

    # Main options
    main_group = parser.add_argument_group('Main options')
    main_group.add_argument("-d", "--db_dir", dest="databasefolder", required=True,
                            help="directory containing MIrROR database [required]", metavar="DBDIR")
    main_group.add_argument("-o", "--output_dir", dest="outputfolder", default='./Result',
                            help="specify directory to output files (default: ./Result)", metavar="OUTDIR")
    main_group.add_argument("-t", "--threads", dest="threads", type=int, default=4,
                            help="number of threads (default: 4)", metavar="INT")

    # Mapping options
    mapping_group = parser.add_argument_group('Mapping options')
    mapping_group.add_argument("-x", "--preset", dest="preset", default="map-ont",
                               help="preset options to optimize alignment for different platforms.\n"
                               "map-pb/map-hifi/map-ont - CLR/HiFi/Nanopore (default: map-ont)", metavar="STR")
    mapping_group.add_argument("-M", "--minibatch", dest="minibatch", default="500M",
                               help="number of query bases loaded to memory at once. "
                               "K/M/G suffix accepted. (default: 500M)", metavar="NUM")

    # Threshold options
    threshold_group = parser.add_argument_group('Threshold options')
    threshold_group.add_argument("-m", "--residuematches", dest="residuemaches", type=int, default=2500,
                                 help="minimum number of residue matches (default: 2500)", metavar="INT")
    threshold_group.add_argument("-b", "--blocklength", dest="blocklength", type=int, default=3500,
                                 help="minimum alignment block length (default: 3500)", metavar="INT")
    threshold_group.add_argument("-n", "--Normalization", dest="normalization", action="store_true",
                                 help="normalize counts by 16S-23S rRNA operon copy number")

    # Visualization options
    viz_group = parser.add_argument_group('Visualization options')
    viz_group.add_argument("-K", "--Krona", dest="krona", action="store_true",
                           help="create a Krona plot (requires KronaTools)")
    viz_group.add_argument("-S", "--Stackedplot", dest="stackedplot", action="store_true",
                           help="create stacked bar plots (requires pandas, matplotlib)")
    viz_group.add_argument("-V", "--visualized", dest="visualized", action="store_true",
                           help="perform all visualization tasks (same as -K -S)")

    # help / version
    other_group = parser.add_argument_group('Others')
    other_group.add_argument("-h", "--help", dest="help", action="help",
                             help='show this help message and exit')
    other_group.add_argument("-v", "--version", action="version",
                             version="%s\n%s" % (MYVERSION, MYBUILT))

    options = parser.parse_args()
    args = options.inputfiles

    if len(args) == 0:
        parser.print_help()
        WriteLog("Error", "INPUTFILE required.")
        sys.exit(1)

    if options.databasefolder is None:
        WriteLog("Error", "Need path to 16S-23S rRNA operon database (-d DBDIR)")
        parser.print_help()
        sys.exit(1)

    if not options.databasefolder.startswith("/"):
        options.databasefolder = GETDIR(os.getcwd(), options.databasefolder)
    if not os.path.isdir(options.databasefolder):
        WriteLog("Error", "Database folder not found: %s" % options.databasefolder)
        parser.print_help()
        sys.exit(1)

    if not options.outputfolder.startswith('/'):
        options.outputfolder = GETDIR(os.getcwd(), options.outputfolder)
    if os.path.isdir(os.path.dirname(options.outputfolder)):
        if os.path.isdir(options.outputfolder):
            WriteLog("Error", "Output folder already exists: %s" % options.outputfolder)
            parser.print_help()
            sys.exit(1)
        else:
            os.mkdir(options.outputfolder)
    else:
        WriteLog("Error", "Parent directory not found: %s" % os.path.dirname(options.outputfolder))
        parser.print_help()
        sys.exit(1)
    
    Main(options, args, sys.argv)
