#!/usr/bin/python

import sys, re, os, datetime, platform, glob, math, time

MYVERSION = 'MIrROR tools v0.1'
BUILDDATE = '2021-01-21'
MYBUILT = "built %s" % (datetime.datetime.strptime(BUILDDATE,'%Y-%m-%d').strftime('%b %d %Y'))
DBVERSION = ''
OSSYSTEM = platform.system()
options = ''
args = ''
mycommand = ''
fplog = None
skipmapping = False

fastaextlist = ['.fasta','.fasta.gz','.fa','.fa.gz']
fastqextlist = ['.fastq','.fastq.gz','.fq','.fq.gz']
pafextlist = ['.paf']
ProcessStep = {1:'filecheck',2:'mapping',3:'classification',4:'OTUtable',5:'graph'}

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
    fploglist = glob.glob("%s/%s*.log" % (options.outputfolder, os.path.basename(options.outputfolder.rstrip("/"))) )
    logindex = []
    for fplogfile in fploglist:
        index = os.path.basename(fplogfile).replace(os.path.basename(options.outputfolder.rstrip("/")),"").replace(".log","")
        if( index.isdigit() ):
            logindex.append(int(index))

    currentindex = 0
    if( len(logindex) > 0 ):
        currentindex = max(logindex)+1

    fplogfile = "%s/%s.log" % (options.outputfolder,os.path.basename(options.outputfolder.rstrip("/")))

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
def ConfirmInput_main():
    rv = True
    inputformat = ''
    inputformatstr = ''
    groupexist = 'No'

    mmiversion = []
    mmilist = glob.glob("%s/MIrROR_DB_*.mmi" % (options.databasefolder))
    for f in mmilist:
        vnum = os.path.basename(f).replace("MIrROR_DB_","").replace(".mmi","")
        mmiversion.append(vnum)
    tsvversion = []
    tsvlist = glob.glob("%s/MIrROR_DB_*.tsv" % (options.databasefolder))
    for f in tsvlist:
        vnum = os.path.basename(f).replace("MIrROR_DB_","").replace(".tsv","")
        tsvversion.append(vnum)
    
    comversion = list(set(mmiversion).intersection(set(tsvversion)))
    comversion.sort(reverse=True)
    if( len(comversion) == 0 ):
        databasestr = "Could not find a suitable file to use as a database."
        rv = False
    else:
        global DBVERSION
        DBVERSION = comversion[0]
        
        database_mmi = "%s/MIrROR_DB_%s.mmi" % (options.databasefolder, DBVERSION)
        database_tsv = "%s/MIrROR_DB_%s.tsv" % (options.databasefolder, DBVERSION)
        
        fp=open(database_tsv,"r")
        headtemp = fp.readline().strip().split()[0].replace("#","")
        dbbuild_date = datetime.datetime.strptime(headtemp+" 00:00:00.00000",'%Y-%m-%d %H:%M:%S.%f')
        fp.close()
                
        databasestr = ''
        databasestr += "%s" % (database_mmi)
        databasestr += "\n%28s %s" % ("",database_tsv)
        databasestr += "\n%28s Database release %s (built %s)" % ("",DBVERSION.replace("r",""),dbbuild_date.strftime('%b %d %Y'))
    
    Readinput()

    if( len(grouporder) > 0 ):
        groupexist = 'Yes'
        
    for sample, infolist in filepersample.items():
        for info in infolist:
            if( inputformatstr == '' ):
                inputformatstr += "[%s] %s" % (info['format'], info['exist'])
            else:
                inputformatstr += "\n%28s [%s] %s" % ("",info['format'], info['exist'])

            if( 'File Not Exist' in info['exist'] ):
                rv = False

    WriteLog("LOG", "%s %s" % (MYVERSION, MYBUILT), True)
    WriteLog("LOG", "%-26s : %s" % ("Command"," ".join(mycommand)), True)
    WriteLog("LOG", "%-26s : %s" % ("DataBase",databasestr), True)
    WriteLog("LOG", "%-26s : %s" %("Input Format",inputformatstr), True)
    WriteLog("LOG", "%-26s : %s" %("Group information",groupexist), True)
    WriteLog("LOG", "%-26s : Number of residue matches-%dbp, Alignment block length - %dbp" %("Read filtering threshold",options.residuemaches, options.blocklength), True)
    WriteLog("LOG", "\n", True)

    return rv

#############
def Readinput():
    global filepersample
    global sampleorder
    global grouporder

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
    outfolder = "%s/%s"%(options.outputfolder, ProcessStep[1])

    savedfile = 0
    rv = True
    for sample, filelist in filepersample.items():
        savedfile += 1
        if( len(filelist) > 1 ):
            WriteLog("Error", "Sample File Duplicate : %s" %(sample))
            rv = False
        else:
            for fileinfo in filelist:
                if( not os.path.isfile(fileinfo['loc']) ):
                    WriteLog("Error", "File not exist : %s" % (fileinfo['loc']))
                    rv = False        
    WriteLog("INFO", "%d samples procedding"%(savedfile))

    return rv

###################
def step2_Mapping(exemapping):
    outfolder = "%s/%s"%(options.outputfolder, ProcessStep[2])

    ref = "%s/MIrROR_DB_%s.mmi" % (options.databasefolder, DBVERSION)
    if( not os.path.isfile(ref)):
            WriteLog("Error","'%s' not exist" %(ref))
            return False

    mustremovefile = []
    systemstr = 'which minimap2 > %s/minimap2.confirm' % (outfolder)

    os.system(systemstr)
    mustremovefile.append("%s/minimap2.confirm"%(outfolder))
    
    fp=open("%s/minimap2.confirm" % (outfolder),"r")
    head = fp.readline()       
    fp.close()
    if( head == '' ):
        WriteLog("Error","'minimap2' not exist in $PATH" )
        for rfile in mustremovefile:
            if( os.path.isfile(rfile)):
                os.unlink(rfile)
        return False

    for rfile in mustremovefile:
        if( os.path.isfile(rfile)):
            os.unlink(rfile)

    mapping_list = []
    idx = 0
    for sample, sampleinfo in filepersample.items():
        if( sampleinfo[0]['format'] in ['FASTQ','FASTA']):
            idx += 1
            WriteLog("INFO","%d / %d Mapping..." % (idx, exemapping))

            mapping_path = "%s/%s_minimap.paf" % (outfolder,os.path.basename(sample))
            log_mapping_path = "%s/%s_minimap.paf.log" % (outfolder,os.path.basename(sample))
            systemstr = "minimap2 --secondary=no -x map-ont -K %s -t %d -c %s %s -o %s 2> %s" % (options.minibatch, options.threads, ref,  sampleinfo[0]['loc'], mapping_path, log_mapping_path)
            
            return_value = os.system(systemstr)
            if( return_value == 0 ):
                mapping_list.append(mapping_path)
                sampleinfo[0]['paf'] = mapping_path

    if( idx != len(mapping_list) ):
        return False

    return True

###################
def step3_Annotaion(skipmapping='False'):
    outfolder = outfolder = "%s/%s/"%(options.outputfolder, ProcessStep[3])

    Annofile = "%s/MIrROR_DB_%s.tsv" % (options.databasefolder, DBVERSION)
    accession_taxa={}
    infile=open(Annofile,"r")

    for full_lineage in infile:
        if( full_lineage.startswith("#")):
            continue
        if( len(full_lineage.strip()) == 0 ):
            continue
        taxa=full_lineage.strip().split("\t")
        
        accession_taxa[taxa[0]]=(taxa[2],taxa[3],taxa[4],taxa[5],taxa[6],taxa[7])
    infile.close()

    classification_list=[]
    idx = 0
    for sample, sampleinfo in filepersample.items():
        idx += 1
        WriteLog("INFO","%d / %d Annotatioin..." % (idx, len(filepersample)))

        mapping_file = sampleinfo[0]['paf']
        paf_file=open(mapping_file,"r")
        read_count={}
        for line in paf_file:
            element=line.split("\t")
            if (int(element[9]) >= int(options.residuemaches)) and (int(element[10]) >= int(options.blocklength)):
                seqname = element[5][:15]
                read_count.setdefault(seqname,0)
                read_count[seqname] += 1
            else: 
                continue
        paf_file.close()

        result={}
        for seqname, count in read_count.items():
            if( seqname in accession_taxa ):
                seqinfo = accession_taxa[seqname]
                result.setdefault(seqinfo,0)
                result[seqinfo] += count

        sorted_result=sorted(result.items(), key=lambda x: x[-1], reverse=True)

        classification_path = "%s/%s.txt" % (outfolder, os.path.basename(sample))
        outfile=open(classification_path,"w")
        for i in sorted_result:
            outfile.write("\t".join(i[0])+"\t"+str(i[1])+"\n")
        outfile.close()
        classification_list.append(classification_path)

    if( len(filepersample) != len(classification_list) ):
        return False

    return True

###################
def Get_Operon():
    Annofile = "%s/MIrROR_DB_%s.tsv" % (options.databasefolder, DBVERSION)
    taxa_count={}
    infile=open(Annofile,"r")

    try:
        for full_lineage in infile:
            if( full_lineage.startswith('#Accession')):
                head = full_lineage.strip().split("\t")
                o_idx = head.index('operon')
                continue
            elif( full_lineage.startswith("#") ):
                continue

            taxa=full_lineage.strip().split("\t")
            operon = int(taxa[o_idx])
            s = taxa[-1]
            taxa_count.setdefault(s,[])
            taxa_count[s].append(operon)

    except ValueError:
        WriteLog("Error","Get_Operon")

    finally:
        infile.close()

    taxa_ave = {}
    for key, values in taxa_count.items():
        taxa_ave.setdefault(key,sum(values)/len(values))
            
    return taxa_ave

###################
def READ_GROUP(sampleorder, wtype='type1'):
    fp=open(options.groupfile,"r")
    head = fp.readline().split()
    group_level = len(head)-1
    GROUP = {}
    for line in fp:
        line_temp = line.strip().split()
        GROUP.setdefault(line_temp[0],line_temp[1:])
    fp.close()

    groupline = []
    for level in range(len(group_level)):
        wlist = []
        if( wtype == 'type1'):
            wlist.append('GroupLevel%d'%(level))
        for sample in sampleorder:

            try:
                g = GROUP[sample][level]
            except KeyError:
                g = 'None'

            wlist.append(g)

        if( wtype != 'type1' ):
            wlist.append('GroupLevel%d'%(level))

        groupline.append("\t".join(wlist)+"\n")

    return groupline

###################
def Make_Merged_Sample_type(outfile_list, tagstr, outfilename):
    SAMPLEDATA = {}
    bacteria_kind = [] 
    for outfile in outfile_list:
        samplename = os.path.basename(outfile).replace(tagstr,"")
        SAMPLEDATA.setdefault(samplename,{})
        fp=open(outfile,"r")
        for line in fp:
            line_temp = line.strip().split()
            SAMPLEDATA[samplename].setdefault(line_temp[0],line_temp[1])
            bacteria_kind.append(line_temp[0])
        fp.close()


    bacteria_kind = list(set(bacteria_kind))
    bacteria_kind.sort()

    grouplines = []
    for level in range(len(grouporder)):
        wlist = []
        levelname = grouporder[level]
        wlist.append(levelname)
        for sample in sampleorder:
            g = 'None'
            try:
                
                g = filepersample[sample][0]['group'][level]
            except IndexError:
                g = 'None'

            wlist.append(g)
        grouplines.append("\t".join(wlist)+"\n")

    fpout=open(outfilename,"w")
    fpout.write("\t".join(['#Name']+sampleorder)+"\n")
    if( len(grouporder) > 0 ):
        for groupline in grouplines:
            fpout.write(groupline)
    for kind in bacteria_kind:
        wlist = [kind]
        for sample in sampleorder:
            if( kind in SAMPLEDATA[sample] ):
                wlist.append(SAMPLEDATA[sample][kind])
            else:
                wlist.append('0')
        fpout.write("\t".join(wlist)+"\n")
    fpout.close()

    fpout=open(outfilename+".forgraph","w")
    fpout.write("\t".join(['#Name']+sampleorder)+"\n")
    for kind in bacteria_kind:
        wlist = [kind]
        for sample in sampleorder:
            if( kind in SAMPLEDATA[sample] ):
                wlist.append(SAMPLEDATA[sample][kind])
            else:
                wlist.append('0')
        fpout.write("\t".join(wlist)+"\n")
    fpout.close()

###################
def Make_Merged_Sample_type1(outfile_list, tagstr, outfilename):
    SAMPLEDATA = {}
    bacteria_kind = [] 
    for outfile in outfile_list:
        samplename = os.path.basename(outfile).replace(tagstr,"")
        SAMPLEDATA.setdefault(samplename,{})
        fp=open(outfile,"r")
        for line in fp:
            line_temp = line.strip().split()
            SAMPLEDATA[samplename].setdefault(line_temp[0],line_temp[1])
            bacteria_kind.append(line_temp[0])
        fp.close()

    bacteria_kind = list(set(bacteria_kind))
    bacteria_kind.sort()

    grouplines = []
    for level in range(len(grouporder)):
        wlist = []
        levelname = grouporder[level]
        wlist.append(levelname)
        for sample in sampleorder:
            g = 'None'
            try:
                g = filepersample[sample][0]['group'][level]
            except IndexError:
                g = 'None'

            wlist.append(g)
        grouplines.append("\t".join(wlist)+"\n")

    fpout=open(outfilename,"w")
    fpout.write("\t".join(['#Name']+sampleorder)+"\n")
    if( len(grouporder) > 0 ) :
        for groupline in grouplines:
            fpout.write(groupline)
    for kind in bacteria_kind:
        wlist = [kind]
        for sample in sampleorder:
            if( kind in SAMPLEDATA[sample] ):
                wlist.append(SAMPLEDATA[sample][kind])
            else:
                wlist.append('0')
        fpout.write("\t".join(wlist)+"\n")
    fpout.close()

###################
def Make_Merged_Sample_type2(outfile_list, tagstr, outfilename):
    SAMPLEDATA = {}
    bacteria_kind = [] 
    for outfile in outfile_list:
        samplename = os.path.basename(outfile).replace(tagstr,"")
        SAMPLEDATA.setdefault(samplename,{})
        fp=open(outfile,"r")
        for line in fp:
            line_temp = line.strip().split()
            SAMPLEDATA[samplename].setdefault(line_temp[0],line_temp[1])
            bacteria_kind.append(line_temp[0])
        fp.close()

    bacteria_kind = list(set(bacteria_kind))
    bacteria_kind.sort()

    fpout=open(outfilename,"w")
    fpout.write("\t".join(sampleorder+['phylum;class;order;family;genus;species'])+"\n")

    for kind in bacteria_kind:
        wlist = []
        for sample in sampleorder:
            if( kind in SAMPLEDATA[sample] ):
                wlist.append(SAMPLEDATA[sample][kind])
            else:
                wlist.append('0')
        wlist.append(kind)
        fpout.write("\t".join(wlist)+"\n")
    fpout.close()

###################
def step4_Makeoutfile():
    infolder = "%s/%s"%(options.outputfolder, ProcessStep[3])
    outfolder = "%s/%s"%(options.outputfolder, ProcessStep[4])
    SEP = ';'
    
    taxa_ave = Get_Operon()

    idx = 0
    for sample, sampleinfo in filepersample.items():
        idx += 1

        classification_path = "%s/%s.txt" % (infolder, os.path.basename(sample))
        outfile_levelremove_path="%s/%s_out_removelevel.txt" % (outfolder, os.path.basename(sample) )
        outfile_path="%s/%s_out.txt" % (outfolder, os.path.basename(sample) )
        
        WDATA_remove = {}
        fpout_remove=open(outfile_levelremove_path,"w")
        WDATA = {}
        fpout=open(outfile_path,"w")
        fp=open(classification_path,"r")
        category=['p','c','o','f','g','s']
        for line in fp:
            line_temp = line.strip().split()
            new_line_temp = [] 
            idx = 0
            for value in line_temp[:-1]:
                value = "%s__%s" % (category[idx],value)
                new_line_temp.append(value)
                idx += 1

            taxa = line_temp[-2]

            count = int(line_temp[-1])
            if( options.normalization ):
                print("**nomalization **")
                count = count / taxa_ave[taxa]

            tempstr = SEP.join(new_line_temp)
            WDATA.setdefault(tempstr,0)
            WDATA[tempstr] += count

            for i in range(1,len(new_line_temp)+1):
                tempstr = SEP.join(new_line_temp[:i])
                WDATA_remove.setdefault(tempstr,0)
                WDATA_remove[tempstr] += count
        
        for tempstr, count in WDATA.items():
            fpout.write(" ".join(map(str,[tempstr,count]))+"\n")
            
        for tempstr, count in WDATA_remove.items():
            fpout_remove.write(" ".join(map(str,[tempstr,count]))+"\n")
            
        WriteLog("INFO", "Create OTU file : %s"%(outfile_path))
        fp.close()
        fpout.close()
        fpout_remove.close()


    outfile_levelremove_list = []
    outfile_list = []
    for sample, sampleinfo in filepersample.items():
        outfile_path="%s/%s_out.txt" % (outfolder, os.path.basename(sample) )
        outfile_list.append(outfile_path)

        outfile_levelremove_path="%s/%s_out_removelevel.txt" % (outfolder, os.path.basename(sample) )
        outfile_levelremove_list.append(outfile_levelremove_path)
        

    merged_file = "%s/OUTPUT_std.txt" % (outfolder)
    Make_Merged_Sample_type(outfile_list, '_out.txt', merged_file) 
    WriteLog("INFO","Created file : %s" %(merged_file))

    merged_file = "%s/OUTPUT_mpa.txt" % (outfolder)
    Make_Merged_Sample_type1(outfile_levelremove_list, '_out_removelevel.txt', merged_file)

    merged_file = "%s/OUTPUT_std_type2.txt" % (outfolder)
    Make_Merged_Sample_type2(outfile_list, '_out.txt', merged_file)


    for rfile in outfile_list + outfile_levelremove_list:
        if( os.path.isfile(rfile)):
            os.unlink(rfile)        

    return True

###################
def stacked_plot(result, taxa_level, group, outfile):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')

    SEP = ";"
    lineage = ["phylum","class","order","family","genus","species"]
    cmap = cm.get_cmap("Spectral")
    df = pd.read_csv(result, sep="\t")
    df['#Name'] = df['#Name'].apply(lambda taxonomy: taxonomy.split(SEP)[lineage.index(taxa_level)])
    df = df.groupby('#Name').sum()
    df = df.div(df.sum()).transpose()
    f = plt.figure()
    plt.title("Stacked bar plot - "+taxa_level, color='black')
    df.plot(kind='bar',stacked=True,rot=90,cmap=cmap, ax=f.gca())
    plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
    plt.xlabel("sample")
    plt.savefig(outfile,bbox_inches='tight')

    return

###################
def step5_graph():
    infolder = "%s/%s"%(options.outputfolder, ProcessStep[4])
    outfolder = "%s/%s"%(options.outputfolder, ProcessStep[5])
    mirror_abspath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    removedlevel_merged_file = "%s/OUTPUT_mpa.txt" % (infolder)
    merged_file = "%s/OUTPUT_std.txt" % (infolder)
    merged_file_forgraph = "%s/OUTPUT_std.txt.forgraph" % (infolder)
    merged_type2_file = "%s/OUTPUT_std_type2.txt" % (infolder)

    mustremovefile = []

    if(options.krona is not None or options.visualized is not None ):
        krona_log_file = "%s/krona.log" % (outfolder)
        if( os.path.isdir("%s/krona" % (outfolder))):
            os.system("rm -rf %s/krona"%(outfolder))
        mustremovefile.append(krona_log_file)
        

    plot_sh = "%s/graph.sh" % (outfolder)
    fpout=open(plot_sh,"w")
    fpout.write("cd %s\n" % (outfolder))
    if(options.krona is not None or options.visualized is not None ):
        fpout.write("%s/lib/OTUsamples2krona.v0.2.2.sh %s > %s 2>&1\n" % (mirror_abspath, merged_type2_file, krona_log_file) )

    fpout.close()
    mustremovefile.append(plot_sh)

    return_value = os.system("sh %s" %(plot_sh))

    if(options.stackedplot is not None or options.visualized is not None ):
        for domain in ['phylum','class','order','family','genus','species']:
            stacked_file = "%s/stacked_%s.png" % (outfolder, domain)
            stacked_plot(merged_file_forgraph, domain ,"group", stacked_file)

            WriteLog("INFO","Created Stacked plot Graphs : %s" % (stacked_file))
    
    mustremovefile.append(merged_file_forgraph)    

    for rfile in mustremovefile:
        if( os.path.isfile(rfile)):
            writelog = False
            if( os.path.basename(rfile) == 'krona.log'):
                WriteLog("INFO", "Log from krona")
                writelog = True
            if( os.path.basename(rfile) == 'lefse.log'):
                WriteLog("INFO", "Log from LEfSe")
                writelog = True
                
            if( writelog ):
                fp=open(rfile,"r")
                for line in fp:
                    WriteLog("LOG",line.strip(),True)
                fp.close()
                
            os.unlink(rfile)

    if(options.krona is not None or options.visualized is not None ):
        if( os.path.isdir("%s/krona" % (outfolder)) ):
            WriteLog("INFO","Created Krona Graphs in  %s/krona" % (outfolder))

    return True

###################
def Proceding(step):
    rv = False

    if( step == 1):
        WriteLog("INFO","##%s STEP START"%(ProcessStep[step]) )
        MakeFolder(step,False)
        rv = step1_Rawdata()

    elif( step == 2 ):
        exemapping = 0
        for sample, sampleinfo in filepersample.items():
            if( sampleinfo[0]['format'] in ['FASTQ','FASTA']):
                exemapping += 1
        if( exemapping ):
            WriteLog("INFO","%s STEP START"%(ProcessStep[step]) )
            MakeFolder(step)
            rv = step2_Mapping(exemapping)
        else:
            rv = True

    elif( step == 3 ):
        WriteLog("INFO","%s STEP START"%(ProcessStep[step]) )
        MakeFolder(step)
        rv = step3_Annotaion(skipmapping)

    elif( step == 4 ):
        WriteLog("INFO","%s STEP START"%(ProcessStep[step]) )
        MakeFolder(step)
        rv = step4_Makeoutfile()

    elif( step == 5 ):
        WriteLog("INFO","%s STEP START"%(ProcessStep[step]) )
        MakeFolder(step)
        rv = step5_graph() 

    return rv

###################
def Main(l_options,l_args, l_command):
    global options
    global args
    global fplog
    global mycommand
    global skipmapping 

    start = time.time()
    options = l_options
    args = l_args
    mycommand = l_command
    fplogfile = MakeLogfile()
    fplog = open(fplogfile,"w")
    
    rv = ConfirmInput_main()
    if( rv ):
        WriteLog("INFO",'MIrROR tools Start')
        CurrentProcess = {1:'filecheck',2:'mapping',3:'classification',4:'OTUtable'}

        if( options.krona or options.stackedplot or options.visualized ):
            CurrentProcess.setdefault(5,'graph')

        for step, foldname in CurrentProcess.items():
            rv = Proceding(step)
            if( not rv ):
                WriteLog("Error","Error occurred during step %s"%(foldname))
                sys.exit()
        WriteLog("INFO",'MIrROR tools Done!')
    else:
        WriteLog("Error",'Input File Error')

    end = time.time()
    elapsed_time = datetime.timedelta(seconds=end-start)
    WriteLog("INFO", "\n%-26s : %s" % ("Total time used",elapsed_time), True)
    fplog.close()
    print("\nMIrROIR tools log file created : %s\n" % (fplogfile))
