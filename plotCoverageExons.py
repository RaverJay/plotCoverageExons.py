#!/usr/bin/env python
# -*- coding: utf-8 -*-
# plotCoverageExons.py -- Sebastian Krautwurst

import os
import sys
import getopt
import numpy
import pysam
import pickle
import matplotlib as mpl
mpl.use('Agg') # do not require X window
from matplotlib import pyplot as plt

def usage():
    print('plotCoverageExons.py Usage:')
    print('./plotCoverageExons.py [-p <geneposition>] [-options] <bamfile1> [<bamfile2> ...]')
    print('')
    print('-p <geneposition>          specify genomic position with: <chromosome>:<start>-<end>')
    print('-z <sizefactor1>,<sf2>,... specify normalization size factors for each bamfile in order')
    print('-s <+|->                   specify strand ( + or - ) to be included in plot title')
    print('-e <start>-<end>,<s-e>,... select feature regions to be plotted exclusively')
    print('-l <label1>,<l2>,...       specify labels for the regions')
    print('-c <condition1>,<c2>,...   specify patterns to group bamfiles by filename')
    print('-n <name1>,<n2>,...        specify names for the condition groups')
    print('-u <uname>                 include unmatched bamfiles under name <uname>')
    print('-a                         average samples of the same condition, default off (plot all samples)')
    print('-L                         set logarithmic scale for the y-axis')
    print('-S                         add smoothed curve (rolling mean of 49-mers, twice applied)')
    print('-f                         enable samtools-like filtering of reads (default no filtering)')
    print('-o <file>                  specify output filename and thus output format')
    print('-P <file>                  only write coverage data to pickle file')
    exit(1)

if len(sys.argv)==1:
    usage()

# getopt
#########
try:
    opts, args = getopt.getopt(sys.argv[1:], "p:z:s:e:l:c:n:u:o:P:aLSfbwm")
except getopt.GetoptError as err:
    # print help information and exit:
    print(err) # will print something like "option -a not recognized"
    usage()

pysamtoggle = True
usewigtoggle = False
usebedtoggle = False
geneposstr = ''
normalizetoggle = False
normalizestr = ''
strandtoggle = False
strandstr = ''
exontoggle = False
exonstr = ''
labeltoggle = False
labelstr = ''
conditiontoggle = False
conditionstr = ''
namestoggle = False
namesstr = ''
unmatchedtoggle = False
unmatchedstr = ''
averagetoggle = False
logscaletoggle = False
smoothingtoggle = False
filtertoggle = False
outfilenametoggle = False
outfilename = ''
pickletoggle = False
picklename = ''
manjatoggle = False


for o, a in opts:
    if (o == "-p"):
        geneposstr = a
    elif (o == "-z"):
        normalizetoggle = True
        normalizestr = a
    elif (o == "-s"):
        strandtoggle = True
        strandstr = a
    elif (o == "-e"):
        exontoggle = True
        exonstr = a
    elif (o == "-l"):
        labeltoggle = True
        labelstr = a
    elif (o == "-c"):
        conditiontoggle = True
        conditionstr = a
    elif (o == "-n"):
        namestoggle = True
        namesstr = a
    elif (o == "-u"):
        unmatchedtoggle = True
        unmatchedstr = a
    elif (o == '-o'):
        outfilenametoggle = True
        outfilename = a
    elif (o == '-P'):
        pickletoggle = True
        picklename = a
    elif (o == '-a'):
        averagetoggle = True
    elif (o == '-L'):
        logscaletoggle = True
    elif (o == '-S'):
        smoothingtoggle = True
    elif (o == '-f'):
        filtertoggle = True
    elif (o == '-b'):
        usebedtoggle = True
        pysamtoggle = False
    elif (o == '-w'):
        usewigtoggle = True
        pysamtoggle = False
    elif (o == '-m'):
        manjatoggle = True
    else:
        print("ERROR: Unknown option '"+o+"'.")
        usage()

print('----------')
print('plotCoverageExons.py')

# input checking
#################

# genepos
getfullreflength = False
if geneposstr == '':
    print('> No gene position given, assuming full reference length.')
    getfullreflength = True

else:
    genepos = geneposstr.split(':')
    if (len(genepos) < 2):
        print('ERROR: Invalid gene position.')
        usage()
        
    genechromo = genepos[0]
    generange = genepos[1].split('-')
    if (len(generange) < 2):
        print('ERROR: Invalid gene position.')
        usage()
        
    genestart = generange[0]
    geneend = generange[1]
    if (not genestart.isdigit() or not geneend.isdigit()):
        print('ERROR: Invalid gene position.')
        usage()

    genestarti = int(genestart)
    geneendi = int(geneend)


# strand
if (strandtoggle):
    if (strandstr == '+'):
        strandstr = 'Plus'
    elif (strandstr == '-'):
        strandstr = 'Minus'
    else:
        print('ERROR: Invalid strand: '+strandstr)
        usage()


# conditions
if (conditiontoggle):
    conditionslist = conditionstr.split(',')
    for condition in conditionslist:
        if (condition == ''):
            print("ERROR: Empty condition.")
            usage()
else:
    conditionslist = []

# names for conditions
namesdict = dict()
if (namestoggle):
    if (not conditiontoggle):
        print('ERROR: no conditions to be named')
        usage()
    nameslist = namesstr.split(',')
    for name in nameslist:
        if (name == ''):
            print("ERROR: Empty condition name.")
            usage()
    if (len(nameslist) != len(conditionslist)):
        print("ERROR: Incorrect amount of condition names.")
        usage()
        
    for i in range(len(conditionslist)):
        namesdict[conditionslist[i]] = nameslist[i]
else:
    for condition in conditionslist:
        namesdict[condition] = condition

# bamfiles
bamfileslist = args

if len(bamfileslist)==0:
    print("ERROR: No bamfiles given!")
    usage()

# expand ~ or ~user
bamfileslist = [os.path.expanduser(path) for path in bamfileslist]

for bamfile in bamfileslist:
    if (not os.path.isfile(bamfile)):
        print("ERROR: '"+bamfile+"' is not a file.")
        usage()

# get reference length
if getfullreflength:
    
    print('> Trying to get the reference length from first bamfile ...')

    pysambam = pysam.AlignmentFile(bamfileslist[0], "rb")
    header = pysambam.header
    #print(header)
    ref = header['SQ'][0]

    genechromo = ref['SN']
    genestarti = 1
    genestart = '1'
    geneendi = ref['LN']
    geneend = str(geneendi)

    if geneendi > 50000:
        # if input('>>>!<<< Reference is longer than 50kb, continue? [y/n]: ') not in ['y', 'Y', 'yes', 'Yes']:
        #     print('> Aborting.')
        #     exit(1)
        print('>>>!<<< Reference is longer than 50kb - this may take a while.')

    geneposstr = genechromo + ':' + genestart + '-' + geneend

# exons
if (exontoggle):
    exonlist = []
    exonstrspl = exonstr.split(',')
    for exon in exonstrspl:
        exonspl = exon.split('-')
        if (len(exonspl) != 2):
            print('ERROR: invalid feature position: '+exon)
            usage()
        if (not exonspl[0].isdigit() or not exonspl[1].isdigit()):
            print('ERROR: invalid feature position: '+exon)
            usage()
        estart = int(exonspl[0])
        eend = int(exonspl[1])
        if (estart > eend):
            print('ERROR: invalid feature position: '+exon)
            usage()
        if (estart < genestarti or estart > geneendi or eend > geneendi or eend < genestarti):
            print('ERROR: feature position not inside gene position.')
            usage()
        if (exonlist != []):
            lasteend = exonlist[-1][1]
            if (estart <= lasteend):
                print('ERROR: features overlapping: '+str(exonlist[-1][0])+'-'+str(exonlist[-1][1])+' and '+exon)
                usage()
        exonlist.append((estart, eend))

    #print(exonlist)
    exonlistrel = [(start-genestarti, end-genestarti) for (start, end) in exonlist]
    #print(exonlistrel)

# labels
if (labeltoggle):
    if (not exontoggle):
        print('ERROR: no features to be labeled')
        usage()
    labellist = []
    labelstrspl = labelstr.split(',')
    if (len(labelstrspl) != len(exonlist)):
        print('ERROR: invalid number of feature labels')
        usage()
    for label in labelstrspl:
        if (label == ''):
            print('ERROR: empty feature label')
            usage()
        labellist.append(label)

print('> Input checks done.')
print('> Gene position: '+geneposstr)

genelengthint = geneendi - genestarti
print('> Length of gene: '+str(genelengthint+1))
if (exontoggle):
    print('> Relative exon regions selected: ', end='')
    for (es,ee) in exonlistrel:
        print(str(es)+'-'+str(ee)+'    ', end='')
    print('')

print('> '+str(len(bamfileslist))+' Bamfile(s)')
#c=0
#for f in bamfileslist:
#    c+=1
#    print('> '+str(c)+': '+f)

#print('> Making coverage files')


# normalization

if (normalizetoggle):
    normalizedict = dict()
    
    normalizelist = normalizestr.split(',')
    if (len(normalizelist) != len(bamfileslist)):
        print("ERROR: Incorrect amount of normalization size factors.")
        usage()
    normalizelistf = []
    for sfactor in normalizelist:
        try:
            sff = float(sfactor)
        except ValueError:
            print("ERROR: Invalid normalization size factor: "+sfactor)
            usage()
        normalizelistf.append(sff)
        
    print('> '+str(len(normalizelistf))+" Normalization size factors: ", end='')
    for f in normalizelist:
        print(f+' ', end='')
    print('')


# misc

outbeddir = '_'.join(['plotCoverage', genechromo, genestart, geneend, 'outbed'])
if (usebedtoggle):
    os.system('mkdir -p '+outbeddir)

outwigdir = '_'.join(['plotCoverage', genechromo, genestart, geneend, 'outwig'])
if (usewigtoggle):
    os.system('mkdir -p '+outwigdir)

covfilelist = []
coveragedata = dict()
sampledict = dict()

# pysam
########

if (pysamtoggle):
    for bamfile in bamfileslist:
        
        baifile = bamfile + '.bai'
        if (not os.path.isfile(baifile)):
            print("ERROR: Missing bam index for file: "+bamfile)
            usage()
        
        pysambam = pysam.AlignmentFile(bamfile, "rb")
        if (filtertoggle):
            columns = pysambam.pileup(genechromo, genestarti-1, geneendi, max_depth=9999999, truncate=True, min_base_quality=0, stepper='all') # 0-based, inclusive, does standard samtools depth like filtering of reads
        else:
            columns = pysambam.pileup(genechromo, genestarti-1, geneendi, max_depth=9999999, truncate=True, min_base_quality=0, stepper='nofilter') # 0-based, inclusive, no filter
                
        data = list()
        runthrough = False
        
        try:
            col = next(columns)
        except StopIteration:
            #print('NO COLUMNS')
            runthrough = True
        
        for posidx in range(genestarti-1, geneendi+1): # include last one
            
            if (runthrough):
                data.append(0)
                continue
            
            while (col.pos < posidx):
                try:
                    col = next(columns)
                except StopIteration:
                    #print('RUNTHROUGH! '+str(posidx))
                    runthrough = True
                    break
                    
            if (runthrough):
                data.append(0)
                continue
                
            if(col.pos > posidx):
                data.append(0)
                continue
                
            elif(col.pos == posidx):
                cnt = 0
                for read in col.pileups:
                    if (read.is_del == 0 and read.is_refskip == 0):
                        cnt += 1
                #print('Pos: '+str(col.pos)+'\tReads: '+str(col.n)+'\tCount: '+str(cnt))
                data.append(cnt)
                
            else:
                # should never happen
                print('SUPERFAIL!')
                exit(4)
                
        #print(data)
        #print(len(data))
        
        
        samplename = bamfile.rsplit('/',1)[-1].rsplit('.',1)[0]
        
        #print(samplename)
        #print(numpy.mean(data))
        coveragedata[samplename] = data
        covfilelist.append(samplename)

# make coverage data
#####################

bam2wig = '/home/ya86gul/programs/bin/bam2wig'
if (usewigtoggle):
    for bamfile in bamfileslist:

        # call bam2wig: -F 0 to disable filtering of secondary alignments etc
        bam2wigcall = bam2wig+' -F 0 -c ' + geneposstr + ' -D ' + outwigdir + ' ' + bamfile + ' 2>/dev/null'
        #print('>>> '+bam2wigcall)
        
        samplename = bamfile.rsplit('/',1)[-1].rsplit('.',1)[0]
        outwigfile = outwigdir + '/' + samplename + '-c' + genechromo + '_' + genestart + '-' + geneend + '-F0.wig'
        covfilelist.append(outwigfile)
        
        os.system(bam2wigcall)

if (usebedtoggle):
    for bamfile in bamfileslist:

        samplefile = bamfile.rsplit('/',1)[-1]
        outcovfile = outbeddir + '/' + samplefile + '.bed'

        covfilelist.append(outcovfile)
        
        samtoolscall = 'samtools depth -a -r ' + geneposstr + ' ' + bamfile + ' >' + outcovfile
        #print('>>> '+samtoolscall)
        os.system(samtoolscall)

# read coverage data
#####################

#print('> Reading coverage data')

if (usewigtoggle):
    for covfile in covfilelist:
        #print(covfile)
        sampledict[covfile] = covfile.rsplit('/',1)[-1].rsplit('.',1)[0]
        with open(covfile,'r') as wigfh:
            coveragedata[covfile] = [0.0 for i in range(genelengthint)]
            lines = wigfh.readlines()
            for line in lines:
                if (line[:12] == 'variableStep'):
                    spanstr = line.strip().split('=')[-1]
                    if (spanstr == ''):
                        break
                    span = int(spanstr)
                else:
                    lsplit = line.strip().split(' ')
                    pos = int(lsplit[0])-genestarti
                    if (pos < 0):
                        span += pos
                        if (span < 0):
                            span = 0
                        pos = 0
                    
                    cov = float(lsplit[1])
                    coveragedata[covfile][pos:pos+span] = [cov for i in range(span)]
                    
                    #print(str(pos)+' to '+str(pos+span)+' is '+str(cov))
        coveragedata[covfile] = coveragedata[covfile][:genelengthint]


if (usebedtoggle):
    for covfile in covfilelist:
        sampledict[covfile] = covfile.rsplit('/',1)[-1].rsplit('.',1)[0]
        with open(covfile,'r') as bedfh:
            coveragedata[covfile] = []
            lines = bedfh.readlines()
            for line in lines:
                coveragedata[covfile].append(float(line.strip().split('\t')[2]))


# pickle output
################

if (pickletoggle):
    print('> Pickle output mode: Writing pickled coverage data to '+picklename)
    with open(picklename, 'wb') as pfh:
        pickle.dump(coveragedata, pfh)
    print('> Done, exiting.')
    exit(0)


# normalization
################

if (normalizetoggle):
    for i in range(len(covfilelist)):
        normalizedict[covfilelist[i]] = normalizelistf[i]
    
    for covfile in covfilelist:
        sf = normalizedict[covfile]
        coveragedata[covfile] = [x*sf for x in coveragedata[covfile]]


# conditions and colors
########################
colors = ['b-','r-','g-','y-','m-','c-','b--','r--','g--','y--','m--','c--']
colordict = dict()

groups = dict()
groupcolordict = dict()
if (conditiontoggle):
    
    for condition in conditionslist:
        groups[condition] = []
        groupcolordict[condition] = colors[0]
        colors = colors[1:] + colors[:1] #rotate
    groupcolordict['UNMATCHED'] = colors[0]

    for covfile in covfilelist:
        for condition in conditionslist:
            if (condition in covfile): # string matching
                groups[condition].append(covfile)
                colordict[covfile] = groupcolordict[condition]
                break
        else: # when for did not break
            if (not 'UNMATCHED' in groups):
                groups['UNMATCHED'] = []
            groups['UNMATCHED'].append(covfile)
            colordict[covfile] = groupcolordict['UNMATCHED']
else:
    for covfile in covfilelist:
        colordict[covfile] = colors[0]
        colors = colors[1:] + colors[:1] #rotate

for condition in groups:
    if (len(groups[condition]) <= 0):
        print("ERROR: No samples in group: "+condition)
        usage()



# averaging
############
from operator import add
averagenumsamples = dict()
if(averagetoggle):
    sdcov = dict()
    averagedata = dict()
    for condition in groups:
        numsamples = len(groups[condition])
        averagenumsamples[condition] = numsamples
        meancovs = []
        for covfile in groups[condition]:
            meancovs.append(numpy.mean(coveragedata[covfile]))
            if (not condition in averagedata):
                averagedata[condition] = coveragedata[covfile]
            else:
                averagedata[condition] = list(map(add, averagedata[condition], coveragedata[covfile]))
        averagedata[condition] = [x/numsamples for x in averagedata[condition]]
        sdcov[condition] = numpy.std(meancovs)
        print('> Samples in condition '+condition+': '+str(averagenumsamples[condition]))


# include unmatched
####################
if (unmatchedtoggle):
    conditionslist.append('UNMATCHED')
    namesdict['UNMATCHED'] = unmatchedstr


# plot coverage
################
outname = ''
if (outfilenametoggle):
    outname = outfilename
if (not outfilenametoggle or os.path.isdir(outname)):
    outname += '_'.join(['plotCoverage', genechromo, genestart, geneend])
    if (exontoggle):
        outname += '_'.join(['_exons']+[ str(s)+'-'+str(e) for (s, e) in exonlistrel])
    if (usewigtoggle):
        outname += '_wig'
    if (usebedtoggle):
        outname += '_bed'
    if (averagetoggle):
        outname += '_average'
    if (filtertoggle):
        outname += '_filter'
    outname += '.pdf'


plt.style.use('fivethirtyeight')
mpl.rcParams['lines.linewidth'] = 2

# bist du mandja???
if (manjatoggle):
    mpl.rcParams['lines.linewidth'] = 15

fig, ax = plt.subplots(figsize=(20, 8))

# slice by exons
if (exontoggle):
    maxy = 0.0
    if (averagetoggle):
        for condition in averagedata:
            tmpdata = []
            for (es, ee) in exonlistrel:
                tmpdata += averagedata[condition][es:ee+1]
            #print(len(averagedata[condition]))
            #print(len(tmpdata))
            averagedata[condition] = tmpdata
            maxy = max(maxy, max(tmpdata))
    else:
        for covfile in covfilelist:
            tmpdata = []
            for (es, ee) in exonlistrel:
                tmpdata += coveragedata[covfile][es:ee+1]
            #print(len(coveragedata[covfile]))
            #print(len(tmpdata))
            coveragedata[covfile] = tmpdata
            maxy = max(maxy, max(tmpdata))



if (exontoggle):
    # vertical bounding lines
    bounds = []
    for (s, e) in exonlistrel:
        if(bounds == []):
            bounds.append(e-s)
        else:
            bounds.append(bounds[-1]+e-s+1)
    #print(bounds)
    #print(len(x))
        
    for b in bounds[:-1]:
        ax.plot((b,b),(0,maxy),color='black')
    
    # x labels
    ticks = [0]
    labels = ['1']
    ecnt = 0
    for b in bounds:
        # midpoint
        ticks.append((ticks[-1]+b)/2)
        ecnt += 1
        if (labeltoggle):
            labels.append('\n'+labellist[ecnt-1])
        else:
            labels.append('\nExon '+str(ecnt))
        # bound
        ticks.append(b)
        labels.append(str(b+1))
    #print(ticks)
    #print(labels)
    plt.xticks(ticks, labels)
    ax.xaxis.grid(False)


# smoothing code
def running_mean(x, N):
    cumsum = numpy.cumsum(numpy.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def edge_sums(values, N):
    pad = abs(N) // 2
    edge = []
    if N > 0:
        for n in range(pad+1, 2*pad+1):
            edge.append(numpy.mean(values[:n]))
    elif N < 0:
        for n in range(-2*pad, -pad):
            edge.append(numpy.mean(values[n:]))
    else:
        print('Bugged smoothing code!')
        exit
    return edge

def smooth_signal(values, N, keepZeroes=False):
    smoo = running_mean(values, N)
    smoo = numpy.insert(smoo, 0, edge_sums(values, N))
    smoo = numpy.append(smoo, edge_sums(values, -N))
    if keepZeroes:
        for i, v in enumerate(values):
            if v == 0.0:
                smoo[i] = 0.0
    return smoo

# smoothing window size
smooN = 49


# plot it!
meancov = dict()
labelwarn = False
if (averagetoggle):
    x = range(len(averagedata[conditionslist[0]]))
    for condition in conditionslist:
        meancov[condition] = numpy.mean(averagedata[condition])
        curlabel = namesdict[condition]+' - #samples: '+str(averagenumsamples[condition])+' - aggregate mean coverage: '+"{:.2f}".format(meancov[condition])+' - sd between samples: '+"{:.2f}".format(sdcov[condition])

        ax.plot(x, averagedata[condition], groupcolordict[condition], label=curlabel, alpha=0.75)
        if smoothingtoggle:
            ax.plot(x, smooth_signal(smooth_signal(averagedata[condition], smooN), smooN))

        ax.fill_between(x, [0.0 for x in range(len(x))], averagedata[condition], color=groupcolordict[condition][0], alpha=0.2)
        if curlabel[0]=='_':
                labelwarn = True
else:
    x = range(len(coveragedata[covfilelist[0]]))
    if (conditiontoggle):
        for condition in conditionslist:
            for covfile in groups[condition]:
                sampledict[covfile] = covfile
                meancov[covfile] = numpy.mean(coveragedata[covfile])

                ax.plot(x, coveragedata[covfile], colordict[covfile], label=sampledict[covfile]+' - mean coverage: '+"{:.2f}".format(meancov[covfile]) , alpha=0.7)
                if smoothingtoggle:
                    ax.plot(x, smooth_signal(smooth_signal(coveragedata[covfile], smooN), smooN))
                    

    else:
        for covfile in covfilelist:
            sampledict[covfile] = covfile
            meancov[covfile] = numpy.mean(coveragedata[covfile])

            ax.plot(x, coveragedata[covfile], colordict[covfile], label=sampledict[covfile]+' - mean coverage: '+"{:.2f}".format(meancov[covfile]) , alpha=0.7)
            if smoothingtoggle:
                ax.plot(x, smooth_signal(smooth_signal(coveragedata[covfile], smooN), smooN))


# log scale on y axis
if (logscaletoggle):
    ax.set_yscale('log')


# set ylims
if (logscaletoggle):
    plt.ylim(0.5, plt.ylim()[1])
else:
    yup = plt.ylim()[1]
    plt.ylim(-yup/40, yup+(yup/40))


# legend
ax.legend(bbox_to_anchor=(1.0, 0.97), loc=4, borderaxespad=0.)
if labelwarn:
    print('> WARNING: Condition names starting with _ will not show up on the legend. Use -n to assign meaningful names to conditions.')




titlestr = 'Coverage for '+geneposstr
if(exontoggle):
    titlestr += ' - specific regions'
if(normalizetoggle):
    titlestr += ' - normalized'
if (strandtoggle):
    titlestr += '\nStrand: '+strandstr
plt.title(titlestr, loc='left')

plt.savefig(outname, bbox_inches='tight')

print('> Figure written to file: '+outname)
print('> All done.')
print('----------')
exit(0)
