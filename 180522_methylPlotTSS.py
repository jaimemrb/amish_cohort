#!/usr/bin/python

'''
The MIT License (MIT)

Copyright (c) 2015 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
]
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''

#pythonTemplate.py <- change to title of your script
#130801 <- date
#Name 


#Description:

#This is a generic python template that has functions from utils.py imported and can be used on CFCE1



#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================

import sys

print "Using python version %s" % sys.version


#importing utils package
sys.path.append('/storage/goodell/home/jmreyes/pipeline_jr/')
import utils
import os
import datetime
import time

#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section


dataFile ='/location/file.txt'
genome = 'hg19'
projectFolder = '/storage/goodell/home/jmreyes/projects/'

#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================


annotFile = '/storage/goodell/home/jmreyes/pipeline/annotation/hg19_refseq.ucsc'

#================================================================================
#===============================MAIN RUN=========================================
#================================================================================

def main():

    #get WGBS files
    projectFolder = '/storage/goodell/projects/jmreyes/amish_ayala/'
    wgbsList = ['dc5_mutant_BS.txt', 'dc3_mutant_BS.txt', 'dc15_WT_BS.txt', 'dc16_WT_BS.txt']

    #import out genes of interest
    sigGenesFile = utils.parseTable(projectFolder+'tables/Amish_significant.txt', '\r')
    sigTable = [x.split('\t') for x in sigGenesFile[0]] 
    sigGenes = [x[0] for x in sigTable]

    #make start dict containing all TSS start sites
    startDict = utils.makeStartDict(annotFile)
    
    #converter form refseq to gene name
    revDict = {}
    for name in startDict.keys():
        revDict[startDict[name]['name']] = name

    #get out subset of genes
    sigLoci = []
    window = 500
    for gene in sigGenes:
        if gene in revDict.keys():
            refSeq = revDict[gene]
            geneChr = startDict[refSeq]['chr']
            geneStart = startDict[refSeq]['start']
            geneEnd = startDict[refSeq]['end']
            geneSense = startDict[refSeq]['sense']
            newLocus = [geneChr, geneStart[0]-window, geneStart[0]+window, geneSense, gene+':'+refSeq]
            sigLoci.append(newLocus)
            
        else:
            refSeq = 'NA'

#    print len(sigLoci)
#    print sigLoci[1:5]
    # utils.unParseTable(sigLoci, projectFolder+'bed/Amish_sigTSS_-500_+500.bed', '\t')
    sortedBed = projectFolder+'bed/Amish_sigTSS_-500_+500.sorted.bed'
    

    binNumber = 200
    ts = time.time()
    timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%Hh%Mm%Ss')

    dateFolder = projectFolder+'scripts/'+datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d')+'/'
    utils.formatFolder(dateFolder, True)
    bedList = [sortedBed]

    for wgbsCalls in wgbsList:
        methylPlotBash = [['#!/usr/bin/bash']]
        catBash = []
        wgbsName = wgbsCalls.split('.')[0]
        outDir = projectFolder+'temp/'+wgbsName + '/'
        outBed = projectFolder+'temp/'+wgbsName + '/bed/'

        utils.formatFolder(outDir, True)
        utils.formatFolder(outBed, True)

        ticker = 0
        for bed in bedList:
            bedName = bed.split('.bed')[0].split('/')[-1]
            splitCmd = 'split -l 1000 %s %s' % (bed, outBed+bedName)
            os.system(splitCmd)

            bedSplitList = [x for x in os.listdir(outBed) if bedName in x]
            catBedList = []

            for bed in bedSplitList:
                ticker += 1
                outName = wgbsName+'_'+bed
                if ticker % 10 == 0:
                    sepMark = '&'
                else:
                    sepMark = '&'
                methylCall = 'python /storage/goodell/home/jmreyes/xwing/methylPlot.py -i %s -b %s -o %s -n %s %s' % (projectFolder+'wgbs/'+wgbsCalls, outBed+bed, outDir+outName, binNumber, sepMark)
                methylPlotBash.append([methylCall])
                
                catBedList.append(outDir+outName)

            catBedListSort = sorted(catBedList)
            catOut = projectFolder+'mapped/'+wgbsName+'_'+bedName+'_'+timestamp+'_avgMethyl.txt'
            catCmd = '#cat %s > %s' % (' '.join(catBedListSort), catOut)
            catBash.append([catCmd])
            rmCmd = ['#rm -rf %s' % (outDir)]
            

        outputBash = methylPlotBash + catBash
            
        utils.unParseTable(outputBash, dateFolder+wgbsName+'_TSS_mapping_'+timestamp+'.sh', '\t')

        
main()

