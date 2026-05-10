#!/usr/bin/env python3

import sys, re, os, datetime, platform, glob, math, time, subprocess

MYVERSION = 'MIrROR v1.1'
BUILDDATE = '2026-05-10'
MYBUILT = "Built on %s" % (datetime.datetime.strptime(BUILDDATE,'%Y-%m-%d').strftime('%b %d %Y'))
DB_RELEASE = '' #e.g. 'r01', 'r02' - derived from MIrROR_DB_rXX.tsv filename
OSSYSTEM = platform.system()
options = ''
args = ''
mycommand = ''
fplog = None
skipmapping = False

fastaextlist = ['.fasta','.fasta.gz','.fa','.fa.gz']
fastqextlist = ['.fastq','.fastq.gz','.fq','.fq.gz']
pafextlist = ['.paf']
ProcessStep = {1:'InputQC',2:'ReadMapping',3:'Classification',4:'FeatureTable',5:'Visualization'}

filepersample = {}
sampleorder = []
grouporder = []

#############
def GETDIR(currentloc,outputfolder):
    temp = outputfolder.strip().split("/")
    for v in temp:
        if( v == '.' or v == '' ):
            pass
        elif( v == '..' ):
            currentloc = os.path.dirname(currentloc)
        else:
            currentloc = "%s/%s" %(currentloc,v)

    return currentloc

#############
def WriteLog(msgtype,msg,onlystring=False):
    currenttime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    writestr = "[%s] %s : %s" % (currenttime,msgtype,msg) 
    if(onlystring):
        writestr = "%s" % (msg)

    print(writestr)
    if( fplog is not None ):
        fplog.write("%s\n" % (writestr))

#############
def MakeFolder(curstep, makefolder=True):
    folder = "%s/%s"%(options.outputfolder, ProcessStep[curstep])
    if( makefolder and not os.path.isdir(folder) ):
        os.mkdir(folder)

    return True

#############
def MakeLogfile():
    fplogfile = "%s/%s.log" % (
                options.outputfolder,
                os.path.basename(options.outputfolder.rstrip("/"))
    )

    return fplogfile

#############
def GET_INPUTFORMAT(fname,dir='.'):
    if(fname.startswith("/") or fname.startswith("..") or fname.startswith(".")):
        dir = ''
    else:
        dir = "%s/" % (dir)

    stat = 'File Not Exist : %s%s' % (dir,fname)
    if( os.path.isfile("%s%s" % (dir,fname)) ):
        stat = '%s%s' % ((dir,fname))
    basename = os.path.basename(fname).lower() 
    if( basename.endswith(tuple(fastaextlist)) ):
        return 'FASTA', stat

    elif( basename.endswith(tuple(fastqextlist)) ):
        return 'FASTQ', stat

    elif( basename.endswith(tuple(pafextlist)) ):
        return 'PAF', stat

    else:
        return 'NONE', stat

#############
def _release_num():
    try:
        return int(DB_RELEASE.lstrip('r'))
    except ValueError:
        return 0

#############
def get_taxonomy_categories():
    if _release_num() >= 2:
        return ['d', 'p', 'c', 'o', 'f', 'g', 's']
    else:
        return ['p', 'c', 'o', 'f', 'g', 's']

#############
def get_lineage_levels():
    if _release_num() >= 2:
        return ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    else:
        return ['phylum', 'class', 'order', 'family', 'genus', 'species']

#############
def db_path(ext):
    new_name = "%s/MIrROR_%s%s"    % (options.databasefolder, DB_RELEASE, ext)
    old_name = "%s/MIrROR_DB_%s%s" % (options.databasefolder, DB_RELEASE, ext)
    if os.path.isfile(new_name):
        return new_name
    return old_name

#############
def ConfirmInput_main():
    rv = True
    inputformatstr = ''
    groupexist = 'No'
    
    DB_BUILD_DATES = {'r02': 'Jan 01 2025'}
    
    tsvlist = glob.glob("%s/MIrROR_DB_r*.tsv" % options.databasefolder)+glob.glob("%s/MIrROR_r*.tsv" % options.databasefolder)
    tsvlist.sort(reverse=True)
    
    if len(tsvlist) == 0:
        databasestr = "Could not find MIrROR_DB_rXX.tsv or MIrROR_rXX.tsv in: %s" % options.databasefolder
        rv = False
    else:
        global DB_RELEASE
        basename_tsv = os.path.basename(tsvlist[0])
        DB_RELEASE = basename_tsv.replace("MIrROR_DB_", "").replace("MIrROR_","").replace(".tsv", "")
        release_label = "release %02d" % _release_num()

        database_tsv = tsvlist[0]
        database_mmi = "%s/MIrROR_%s.mmi" % (options.databasefolder, DB_RELEASE)

        if not os.path.isfile(database_mmi):
            database_mmi = db_path(".mmi")

        if DB_RELEASE in DB_BUILD_DATES:
            built_str = DB_BUILD_DATES[DB_RELEASE]
        else:
            try:
                with open(database_tsv, "r") as fp:
                    headtemp = fp.readline().strip().split()[0].replace("#", "")
                dbbuild_date = datetime.datetime.strptime(
                        headtemp + " 00:00:00.00000", '%Y-%m-%d %H:%M:%S.%f')
                built_str = dbbuild_date.strftime('%b %d %Y')
            except (ValueError, IndexError):
                built_str = "unknown"

        if os.path.isfile(database_mmi):
            databasestr = "%s" % database_mmi
            databasestr += "\n%28s %s" % ("", database_tsv)
        else:
            databasestr = "(no .mmi found — mapping will be skipped or use existing PAF)"
            databasestr += "\n%28s %s" % ("", database_tsv)

        databasestr += "\n%28s Database %s (built %s)" % ("", release_label, built_str)

    Readinput()

    if len(grouporder) > 0:
        groupexist = 'Yes'

    for sample, infolist in filepersample.items():
        for info in infolist:
            entry = "[%s] %s" % (info['format'], info['exist'])
            if inputformatstr == '':
                inputformatstr += entry
            else:
                inputformatstr += "\n%28s %s" % ("", entry)
            if 'File Not Exist' in info['exist']:
                rv = False

    WriteLog("LOG", "%s"  % (MYVERSION), True)
    WriteLog("LOG", "%-26s : %s" % ("Command", " ".join(mycommand)), True)
    WriteLog("LOG", "%-26s : %s" % ("DataBase", databasestr), True)
    WriteLog("LOG", "%-26s : %s" % ("Input Format", inputformatstr), True)
    WriteLog("LOG", "%-26s : %s" % ("Group information", groupexist), True)
    WriteLog("LOG", "%-26s : Number of residue matches - %d bp, Alignment block length - %d bp" % (
             "Read filtering threshold", options.residuemaches, options.blocklength), True)
    WriteLog("LOG", "\n", True)

    return rv

#############
def Readinput():
    global filepersample, sampleorder, grouporder

    for argsfile in args:
        basename = os.path.basename(argsfile)
        fullname = argsfile
        if( not fullname.startswith("/")):
            fullname = GETDIR(os.getcwd(),argsfile)
        if( basename.lower().endswith(tuple(fastaextlist+fastqextlist+pafextlist)) ):
            filepersample.setdefault(basename,[])
            format, exist = GET_INPUTFORMAT(fullname,dir='.')
            paf = ''
            if(format == 'PAF'):
                paf = fullname
                
            temp = {}
            temp.setdefault('loc',fullname)
            temp.setdefault('group',[])
            temp.setdefault('format',format)
            temp.setdefault('paf',paf)
            temp.setdefault('exist',exist)
            filepersample[basename].append(temp)

            sampleorder.append(basename)

        else:
            fullname = argsfile
            if( not fullname.startswith("/")):
                fullname = GETDIR(os.getcwd(),argsfile)

            dir = os.path.dirname(fullname)

            fp=open(fullname,"r")
            fileco = 0
            head = fp.readline()
            grouporder = head.strip().split()[1:]
            for line in fp:
                sample = line.strip().split()
                filepersample.setdefault(sample[0],[])
                format, exist = GET_INPUTFORMAT(sample[0],dir)
                paf = ''
                if(format == 'PAF'):
                    paf = "%s/%s"%(dir,sample[0])
                temp = {}
                temp.setdefault('loc',"%s/%s"%(dir,sample[0]))
                temp.setdefault('group',sample[1:])
                temp.setdefault('format',format)
                temp.setdefault('paf',paf)
                temp.setdefault('exist',exist)
                filepersample[sample[0]].append(temp)
                fileco += 1
                sampleorder.append(sample[0])
            fp.close()

    sampleorder.sort()
    return

###################
def step1_Rawdata():
    rv = True
    savedfile = 0
    
    for sample, filelist in filepersample.items():
        savedfile += 1
        if( len(filelist) > 1 ):
            WriteLog("Error", "Sample file duplicated in input: %s" %(sample))
            rv = False
        else:
            for fileinfo in filelist:
                if( not os.path.isfile(fileinfo['loc']) ):
                    WriteLog("Error", "File not found: %s" % (fileinfo['loc']))
                    rv = False

    label = "sample" if savedfile == 1 else "samples"
    WriteLog("INFO", "%d %s to process" % (savedfile, label))

    return rv

###################
def step2_Mapping(exemapping):
    outfolder = "%s/%s"%(options.outputfolder, ProcessStep[2])

    ref = db_path(".mmi")
    if( not os.path.isfile(ref)):
        WriteLog("Error","Database index not found: %s" %(ref))
        return False
    
    try:
        result = subprocess.run(["which", "minimap2"], capture_output=True, text=True, check=True)
        if not result.stdout.strip():
            raise FileNotFoundError
    except (subprocess.CalledProcessError, FileNotFoundError):
        WriteLog("Error", "'minimap2' not found in $PATH. "
                 "Install with: conda install -c bioconda minimap2")
        return False

    mapping_list = []
    idx = 0
    for sample, sampleinfo in filepersample.items():
        if sampleinfo[0]['format'] in ['FASTQ', 'FASTA']:
            idx += 1
            WriteLog("INFO", "%d / %d Mapping..." % (idx, exemapping))

            mapping_path = "%s/%s_minimap.paf" % (outfolder, os.path.basename(sample))
            log_mapping_path = "%s.log" % mapping_path

            cmd = ["minimap2", "--secondary=no", "-x", options.preset,
                   "-K", options.minibatch, "-t", str(options.threads), "-c",
                   ref, sampleinfo[0]['loc'], "-o", mapping_path]

            with open(log_mapping_path, "w") as log_fp:
                ret = subprocess.run(cmd, stderr=log_fp)

            if ret.returncode == 0:
                mapping_list.append(mapping_path)
                sampleinfo[0]['paf'] = mapping_path
            else:
                WriteLog("Error", "Minimap2 failed for sample: %s" % sample)

    return idx == len(mapping_list)

###################
def step3_Annotaion(skipmapping='False'):
    outfolder = "%s/%s/"%(options.outputfolder, ProcessStep[3])

    Annofile = db_path(".tsv")
    accession_taxa={}
    
    with open(Annofile, "r") as infile:
        for full_lineage in infile:
            line = full_lineage.strip()
            if not line or line.startswith("#"):
                continue
            taxa = line.split("\t")
            if taxa[0].lower() in ("accession", "#accession"):
                continue
            if len(taxa) >= 8:
                accession_taxa[taxa[0]] = tuple(taxa[2:])

    classification_list = []
    idx = 0
    for sample, sampleinfo in filepersample.items():
        idx += 1
        WriteLog("INFO", "%d / %d Annotation..." % (idx, len(filepersample)))

        read_count = {}
        with open(sampleinfo[0]['paf'], "r") as paf_file:
            for line in paf_file:
                element = line.split("\t")
                if len(element) < 11:
                    continue
                try:
                    if int(element[9]) >= int(options.residuemaches) and \
                       int(element[10]) >= int(options.blocklength):
                           seqname = element[5][:15]
                           read_count[seqname] = read_count.get(seqname, 0) + 1
                except ValueError:
                    continue

        result = {}
        for seqname, count in read_count.items():
            if seqname in accession_taxa:
                seqinfo = accession_taxa[seqname]
                result[seqinfo] = result.get(seqinfo, 0) + count
        sorted_result = sorted(result.items(), key=lambda x: x[-1], reverse=True)
        classification_path = "%s/%s.txt" % (outfolder, os.path.basename(sample))
        with open(classification_path, "w") as outfile:
            for i in sorted_result:
                outfile.write("\t".join(i[0]) + "\t" + str(i[1]) + "\n")
        classification_list.append(classification_path)

    return len(filepersample) == len(classification_list)

###################
def Get_Operon():
    Annofile = db_path(".tsv")
    taxa_count={}
    o_idx = None
    
    with open(Annofile, "r") as infile:
        header_parsed = False
        o_idx = None
        s_idx = None
        for full_lineage in infile:
            line = full_lineage.strip()
            if not line:
                continue
            if line.startswith('#'):
                cols = [c.replace("#", "").strip().lower() for c in line.split("\t")]
                if 'operon' in cols:
                    o_idx = cols.index('operon')
                    s_idx = len(cols)
                    header_parsed = True
                continue

            cols = line.split("\t")

            if not header_parsed and cols[0].lower() == 'accession':
                col_lower = [c.lower() for c in cols]
                o_idx = col_lower.index('operon')
                s_idx = len(col_lower) - 1
                header_parsed = True
                continue

            try:
                operon = int(cols[o_idx])
                species = cols[s_idx]
                taxa_count.setdefault(species, []).append(operon)
            except (ValueError, IndexError, TypeError):
                continue
    return {key: sum(vals) / len(vals) for key, vals in taxa_count.items()}

###################
def build_grouplines():
    grouplines = []
    for level in range(len(grouporder)):
        wlist = [grouporder[level]]
        for sample in sampleorder:
            g = 'None'
            try:
                g = filepersample[sample][0]['group'][level]
            except IndexError:
                pass
            wlist.append(g)
        grouplines.append("\t".join(wlist) + "\n")

    return grouplines

###################
def load_sample_data(outfile_list, tagstr):
    SAMPLEDATA = {}
    all_taxa = []
    for outfile in outfile_list:
        samplename = os.path.basename(outfile).replace(tagstr, "")
        SAMPLEDATA[samplename] = {}
        with open(outfile, "r") as fp:
            for line in fp:
                parts = line.strip().split()
                if len(parts) >= 2:
                    SAMPLEDATA[samplename][parts[0]] = parts[1]
                    all_taxa.append(parts[0])

    return SAMPLEDATA, sorted(set(all_taxa))

###################
def Make_Merged_Sample_type(outfile_list, tagstr, outfilename):
    SAMPLEDATA, taxa = load_sample_data(outfile_list, tagstr)
    grouplines = build_grouplines()

    with open(outfilename, "w") as fpout:
        fpout.write("\t".join(['#Name'] + sampleorder) + "\n")
        for gl in grouplines:
            fpout.write(gl)
        for kind in taxa:
            row = [kind] + [SAMPLEDATA.get(s, {}).get(kind, '0') for s in sampleorder]
            fpout.write("\t".join(row) + "\n")

    if options.stackedplot or options.visualized:
        with open(outfilename + ".forgraph", "w") as fpout:
            fpout.write("\t".join(['#Name'] + sampleorder) + "\n")
            for kind in taxa:
                row = [kind] + [SAMPLEDATA.get(s, {}).get(kind, '0') for s in sampleorder]
                fpout.write("\t".join(row) + "\n")

###################
def Make_Merged_Sample_type1(outfile_list, tagstr, outfilename):
    SAMPLEDATA, taxa = load_sample_data(outfile_list, tagstr)
    grouplines = build_grouplines()

    with open(outfilename, "w") as fpout:
        fpout.write("\t".join(['#Name'] + sampleorder) + "\n")
        for gl in grouplines:
            fpout.write(gl)
        for kind in taxa:
            row = [kind] + [SAMPLEDATA.get(s, {}).get(kind, '0') for s in sampleorder]
            fpout.write("\t".join(row) + "\n")

###################
def Make_Merged_Sample_type2(outfile_list, tagstr, outfilename):
    SAMPLEDATA, taxa = load_sample_data(outfile_list, tagstr)

    with open(outfilename, "w") as fpout:
        fpout.write("\t".join(sampleorder + ['taxonomy']) + "\n")
        for kind in taxa:
            row = [SAMPLEDATA.get(s, {}).get(kind, '0') for s in sampleorder] + [kind]
            fpout.write("\t".join(row) + "\n")

###################
def step4_Makeoutfile():
    infolder = "%s/%s"%(options.outputfolder, ProcessStep[3])
    outfolder = "%s/%s"%(options.outputfolder, ProcessStep[4])
    SEP = ';'
    
    taxa_ave = Get_Operon()
    category = get_taxonomy_categories()

    outfile_list = []
    outfile_levelremove_list = []

    for sample, sampleinfo in filepersample.items():
        classification_path = "%s/%s.txt" % (infolder, os.path.basename(sample))
        outfile_path = "%s/%s_out.txt" % (outfolder, os.path.basename(sample))
        outfile_levelremove_path = "%s/%s_out_removelevel.txt" % (outfolder, os.path.basename(sample))

        WDATA = {}
        WDATA_remove = {}

        with open(classification_path, "r") as fp:
            for line in fp:
                line_temp = line.strip().split()
                if not line_temp:
                    continue
                taxa_fields = line_temp[:-1]
                try:
                    count = int(float(line_temp[-1]))
                except ValueError:
                    continue

                new_line_temp = [
                        "%s__%s" % (category[i], v) if i < len(category) else v
                        for i, v in enumerate(taxa_fields)
                ]

                taxa_species = line_temp[-2] if len(line_temp) >= 2 else ''
                if options.normalization and taxa_species in taxa_ave and taxa_ave[taxa_species] > 0:
                    count = count / taxa_ave[taxa_species]

                tempstr = SEP.join(new_line_temp)
                WDATA[tempstr] = WDATA.get(tempstr, 0) + count

                for i in range(1, len(new_line_temp) + 1):
                    t2 = SEP.join(new_line_temp[:i])
                    WDATA_remove[t2] = WDATA_remove.get(t2, 0) + count

        with open(outfile_path, "w") as fp:
            for k, v in WDATA.items():
                fp.write("%s %s\n" % (k, round(v)))
        with open(outfile_levelremove_path, "w") as fp:
            for k, v in WDATA_remove.items():
                fp.write("%s %s\n" % (k, round(v)))

        outfile_list.append(outfile_path)
        outfile_levelremove_list.append(outfile_levelremove_path)

    merged_file = "%s/OUTPUT_std.txt" % outfolder
    Make_Merged_Sample_type(outfile_list, '_out.txt', merged_file)
    WriteLog("INFO", "Created feature table (standard version) : %s" % merged_file)

    merged_mpa = "%s/OUTPUT_mpa.txt" % outfolder
    Make_Merged_Sample_type1(outfile_levelremove_list, '_out_removelevel.txt', merged_mpa)
    WriteLog("INFO", "                      (MetaPhlAn version): %s" % merged_mpa)

    merged_type2 = "%s/OUTPUT_std_type2.txt" % outfolder
    Make_Merged_Sample_type2(outfile_list, '_out.txt', merged_type2)
    WriteLog("INFO", "                      (for Krona plot)   : %s" % merged_type2)

    for rfile in outfile_list + outfile_levelremove_list:
        if os.path.isfile(rfile):
            os.unlink(rfile)

    return True

###################
def stacked_plot(result, taxa_level, outfile_pdf):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd

    SEP = ";"
    levels = get_lineage_levels()
    rank_prefixes = [l[0] + "__" for l in levels]

    try:
        cmap = matplotlib.colormaps.get_cmap("Spectral")
    except AttributeError:
        import matplotlib.cm as cm
        cmap = cm.get_cmap("Spectral")

    df = pd.read_csv(result, sep="\t")

    df = df[df['#Name'].str.startswith(tuple(rank_prefixes))]

    target_prefix = taxa_level[0] + "__"

    def extract_level(taxonomy):
        for part in taxonomy.split(SEP):
            if part.startswith(target_prefix):
                return part
        return taxonomy.split(SEP)[0]

    df['#Name'] = df['#Name'].apply(extract_level)
    df = df.groupby('#Name').sum(numeric_only=True)

    col_sums = df.sum()
    col_sums = col_sums.replace(0, 1)
    df = df.div(col_sums).transpose()

    n_samples = len(df)
    n_taxa = len(df.columns)
    fig_width = max(6, n_samples * 0.7 + 3)
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    df.plot(kind='bar', stacked=True, rot=45, cmap=cmap, ax=ax)
    ax.set_title("Stacked bar plot — %s level" % taxa_level, fontsize=12)
    ax.set_xlabel("Sample", fontsize=10)
    ax.set_ylabel("Relative abundance", fontsize=10)
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), ncol=1,
              fontsize=7, title=taxa_level, title_fontsize=8)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
    plt.tight_layout()
    plt.savefig(outfile_pdf, format='pdf', bbox_inches='tight')
    plt.close(fig)

###################
def step5_graph():
    infolder = "%s/%s"%(options.outputfolder, ProcessStep[4])
    outfolder = "%s/%s"%(options.outputfolder, ProcessStep[5])
    mirror_abspath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    merged_file_forgraph = "%s/OUTPUT_std.txt.forgraph" % infolder
    merged_type2_file = "%s/OUTPUT_std_type2.txt" % infolder

    mustremovefile = []

    if options.krona or options.visualized:
        krona_log_file = "%s/krona.log" % outfolder
        krona_dir = "%s/krona" % outfolder
        if os.path.isdir(krona_dir):
            subprocess.run(["rm", "-rf", krona_dir])
        mustremovefile.append(krona_log_file)

        plot_sh = "%s/graph.sh" % outfolder
        with open(plot_sh, "w") as fpout:
            fpout.write("cd %s\n" % outfolder)
            fpout.write("%s/lib/OTUsamples2krona.v0.2.2.sh %s > %s 2>&1\n" % (
                        mirror_abspath, merged_type2_file, krona_log_file))
        mustremovefile.append(plot_sh)

        ret = subprocess.run(["sh", plot_sh])
        if ret.returncode != 0:
            WriteLog("Warning", "Krona plot generation may have encountered errors. "
                     "Ensure KronaTools (ktImportText) is installed and in $PATH.")

    if options.stackedplot or options.visualized:
        levels = get_lineage_levels()
        for level in levels:
            outfile_pdf = "%s/stacked_%s.pdf" % (outfolder, level)
            try:
                stacked_plot(merged_file_forgraph, level, outfile_pdf)
                pad = " " * max(1, 10 - len(level))
                if level == levels[0]:
                    WriteLog("INFO", "Created stacked bar plot (%s)%s: %s" % (level, pad, outfile_pdf))
                else:
                    WriteLog("INFO", "                         (%s)%s: %s" % (level, pad, outfile_pdf))
            except Exeption as e:
                WriteLog("Warning", "Could not generate stacked plot for %s: %s" % (level, str(e)))

    mustremovefile.append(merged_file_forgraph)

    for rfile in mustremovefile:
        if os.path.isfile(rfile):
            if os.path.basename(rfile) == 'krona.log':
                WriteLog("INFO", "Log from OTUsamples2krona")
                with open(rfile, "r") as fp:
                    for line in fp:
                        WriteLog("LOG", line.strip(), True)
            os.unlink(rfile)

    if options.krona or options.visualized:
        if os.path.isdir("%s/krona" % outfolder):
            WriteLog("INFO", "Created Krona plot : %s/krona" % outfolder)

    return True

###################
def Proceeding(step):
    if( step == 1):
        WriteLog("INFO","[STEP 1] %s"%(ProcessStep[step]) )
        MakeFolder(step,False)
        return step1_Rawdata()

    elif( step == 2 ):
        exemapping = sum(
            1 for si in filepersample.values()
            if si[0]['format'] in ['FASTQ', 'FASTA']
        )
        if exemapping:
            WriteLog("INFO", "[STEP 2] %s" % ProcessStep[step])
            MakeFolder(step)
            return step2_Mapping(exemapping)
        else:
            WriteLog("INFO", "[STEP 2] %s skipping..." % ProcessStep[step])
        return True
    
    elif( step == 3 ):
        WriteLog("INFO","[STEP 3] %s"%(ProcessStep[step]) )
        MakeFolder(step)
        return step3_Annotaion(skipmapping)

    elif( step == 4 ):
        WriteLog("INFO","[STEP 4] %s"%(ProcessStep[step]) )
        MakeFolder(step)
        return step4_Makeoutfile()

    elif( step == 5 ):
        WriteLog("INFO","[STEP 5] %s"%(ProcessStep[step]) )
        MakeFolder(step)
        return step5_graph() 

    return False

###################
def Main(l_options,l_args, l_command):
    global options, args, fplog, mycommand, skipmapping

    start = time.time()
    options = l_options
    args = l_args
    mycommand = l_command
    
    fplogfile = MakeLogfile()
    fplog = open(fplogfile,"w")
    
    rv = ConfirmInput_main()
    if( rv ):
        WriteLog("INFO",'MIrROR analysis Start')
        CurrentProcess = {1:'Input QC',2:'Read Mapping',3:'Taxonomic Classification',4:'Feature Table Construction'}

        if( options.krona or options.stackedplot or options.visualized ):
            CurrentProcess[5] = 'Visualization'

        for step, foldname in CurrentProcess.items():
            rv = Proceeding(step)
            if( not rv ):
                WriteLog("Error","Error occurred during step %s"%(foldname))
                fplog.close()
                sys.exit(1)
        WriteLog("INFO",'MIrROR analysis Completed!')
    else:
        WriteLog("Error",'Input File Error -check file paths and database.')

    elapsed = datetime.timedelta(seconds=time.time()-start)
    WriteLog("INFO", "\n%-26s : %s" % ("Total time used",elapsed), True)
    fplog.close()
    print("\nMIrROR analysis log file is created : %s\n" % (fplogfile))
