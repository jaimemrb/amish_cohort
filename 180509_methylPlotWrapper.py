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



def writeSplitBeds(bed, analysisName, outputFolder, window=50, centered = False):

    dmrList = [x for x in utils.parseTable(bed, '\t')]
#    print len(refGenes)

    dmrBed = []
    endsBed = []
    startsBed = []
    centeredBed = []

    if centered == False:

        for line in dmrList:
            dmrID = line[3]

            dmrCoords = [line[0], int(line[1]), int(line[2]), dmrID]
            
            dmrBed.append(dmrCoords)
            startExtend =  [line[0], int(line[1])-window, int(line[1]), dmrID]
            endExtend =  [line[0], int(line[2]), int(line[2])+window,dmrID]
            endsBed.append(endExtend)
            startsBed.append(startExtend)

        print len(dmrBed)
        utils.unParseTable(dmrBed, outputFolder+analysisName+'_BODY_-0_+0.bed', '\t')
        print len(startsBed)
        utils.unParseTable(startsBed, outputFolder+analysisName+'_UPSTREAM_-'+str(window)+'_+'+str(window)+'.bed', '\t')
        print len(endsBed)
        utils.unParseTable(endsBed, outputFolder+analysisName+'_DOWNSTREAM_-'+str(window)+'_+'+str(window)+'.bed', '\t')


    elif centered == True:
        for line in dmrList:
            dmrID = line[3]

            center = (int(line[1])+int(line[2]))/2            
            centeredBed.append([line[0], center - window, center + window, dmrID])

        utils.unParseTable(centeredBed, outputFolder+analysisName+'_CENTERED_-'+str(window)+'_+'+str(window)+'.bed', '\t')        

#    return [outputFolder+'bed/'+analysisName+'_BODY_-0_+0.bed', outputFolder+'bed/'+analysisName+'_UPSTREAM_-'+str(window)+'_+0.bed', outputFolder+'bed/'+analysisName+'_DOWNSTREAM_-0_+'+str(window)+'.bed']



#================================================================================
#===============================MAIN RUN=========================================
#================================================================================

def main():

    #get WGBS files
    projectFolder = '/storage/goodell/projects/jmreyes/amish_ayala/'
    wgbsList = ['dc5_mutant_BS.txt', 'dc3_mutant_BS.txt', 'dc15_WT_BS.txt', 'dc16_WT_BS.txt']
#    bedList = ['hypoDMRsWT.vs.Mut.sorted.bed', 'hyperDMRsWT.vs.Mut.sorted.bed']
#    sortedList = ['TBRS_DMR_hyper.sorted.bed', 'TBRS_DMR_hypo.sorted.bed', 'AML_DMR_hypo.sorted.bed', 'AML_DMR_hyper.sorted.bed']

    binNumber = 200
    ts = time.time()
    timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%Hh%Mm%Ss')

    dateFolder = projectFolder+'scripts/'+datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d')+'/'
    utils.formatFolder(dateFolder, True)

    for bed in sortedList:
        bedPath = projectFolder+'bed/sorted/'+bed
        bedName = bed.split('.bed')[0]+'_extend_'+str(binNumber)

        writeSplitBeds(bedPath, bedName, projectFolder+'bed/sorted/', window=3000, centered = True)
#    keyWords = 
    extendedList = [x for x in os.listdir(projectFolder+'bed/sorted/') if '_DMR_' in x and 'extend' in x and str(binNumber) in x]
    print extendedList

    bedList = extendedList
    methylBashList = []

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
            bedName = bed.split('.bed')[0]
            splitCmd = 'split -l 1000 %s %s' % (projectFolder+'bed/sorted/'+bed, outBed+bedName)
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
            
#        utils.unParseTable(outputBash, dateFolder+wgbsName+'_WTvsMutDMRs_extension_mapping_'+timestamp+'.sh', '\t')
    methylBashList.append([wgbsName+'_mapping.sh'])
        
main()

