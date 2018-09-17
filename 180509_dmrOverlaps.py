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
import random


#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section


dataFile ='/location/file.txt'
genome = 'hg19'
projectFolder = '/storage/goodell/home/jmreyes/projects/amish_ayala/'

#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================

#write your specific functions here
def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here


def main():

    '''
    this is the main run function for the script
    all of the work should occur here, but no functions should be defined here
    '''
    projectFolder = '/storage/goodell/home/jmreyes/projects/amish_ayala/'
    
    #gather up DMR tables
    #ayala MUT vs WT
    mutWT_hypo = utils.parseTable(projectFolder+'bed/hypoDMRsWT.vs.Mut.bed', '\t')
    mutWT_hyper = utils.parseTable(projectFolder+'bed/hyperDMRsWT.vs.Mut.bed', '\t')
    
    mutWT_control = utils.parseTable(projectFolder+'bed/Control_nonDMRsWT.vs.Mut.bed', '\t')



    #ley all
    tbrs_all = utils.parseTable(projectFolder+'bed/TBRS_DMRs.bed', '\t')
    aml_all = utils.parseTable(projectFolder+'bed/AML_DMRs.bed', '\t')

    tbrs_hypo = []
    tbrs_hyper = []

    aml_hypo = []
    aml_hyper = []

    tbrs_all_loci = []
    aml_all_loci = []

    for line in tbrs_all:
        chrom = 'chr'+line[0]
        start = line[1]
        end = line[2]
        tbrs_all_loci.append(utils.Locus(chrom, start, end, '.', 'tbrs_all_'+str(chrom)+':'+str(start)+'-'+str(end)))
        
        if 'hypo' in line:
            chrom = 'chr'+line[0]
            start = line[1]
            end = line[2]
            tbrs_hypo.append(utils.Locus(chrom, start, end, '.', 'tbrs_hypo_'+str(chrom)+':'+str(start)+'-'+str(end)))
        if 'hyper' in line:
            chrom = 'chr'+line[0]
            start = line[1]
            end = line[2]
            tbrs_hyper.append(utils.Locus(chrom, start, end, '.', 'tbrs_hyper_'+str(chrom)+':'+str(start)+'-'+str(end)))

    for line in aml_all:
        chrom = 'chr'+line[0]
        start = line[1]
        end = line[2]
        aml_all_loci.append(utils.Locus(chrom, start, end, '.', 'tbrs_all_'+str(chrom)+':'+str(start)+'-'+str(end)))

        if 'hypo' in line:
            chrom = 'chr'+line[0]
            start = line[1]
            end = line[2]
            aml_hypo.append(utils.Locus(chrom, start, end, '.', 'aml_hypo_'+str(chrom)+':'+str(start)+'-'+str(end)))
        if 'hyper' in line:
            chrom = 'chr'+line[0]
            start = line[1]
            end = line[2]
            aml_hyper.append(utils.Locus(chrom, start, end, '.', 'aml_hyper_'+str(chrom)+':'+str(start)+'-'+str(end)))


    #make locus collection of mutDMRs
    mutWT_control_loci = []

    for line in mutWT_control:

        chrom = line[0]
        start = line[1]
        end = line[2]
        sense = '.'
        locusID = 'control_'+str(chrom)+':'+str(start)+'-'+str(end)
        ratio = line[7]
        if -0.1 < float(ratio) < 0.1:
            new_line = utils.Locus(chrom, start, end, '.', locusID)
            mutWT_control_loci.append(new_line)

    print 'MuitWT control ratio selected:',len(mutWT_control_loci)

    mutWT_hypo_loci = []

    for line in mutWT_hypo:
        chrom = line[0]
        start = line[1]
        end = line[2]
        sense = '.'
        locusID = 'hypo_'+str(chrom)+':'+str(start)+'-'+str(end)
        new_line = utils.Locus(chrom, start, end, '.', locusID)
        mutWT_hypo_loci.append(new_line)

    mutWT_hyper_loci = []

    for line in mutWT_hyper:
        chrom = line[0]
        start = line[1]
        end = line[2]
        sense = '.'
        locusID = 'hyper_'+str(chrom)+':'+str(start)+'-'+str(end)
        new_line = utils.Locus(chrom, start, end, '.', locusID)
        mutWT_hyper_loci.append(new_line)

    #mutWT_hypo_LC = utils.LocusCollection(mutWT_hypo_loci)
    mutWT_all_loci = mutWT_hyper_loci + mutWT_hypo_loci
    
    mutWT_hypo_LC = utils.LocusCollection(mutWT_hypo_loci)

   tbrs_hypo_LC = utils.LocusCollection(tbrs_hypo)
   tbrs_hyper_LC = utils.LocusCollection(tbrs_hyper)
   tbrs_all_LC = utils.LocusCollection(tbrs_all_loci)

   aml_hypo_LC = utils.LocusCollection(aml_hypo)
   aml_hyper_LC = utils.LocusCollection(aml_hyper)
   aml_all_LC = utils.LocusCollection(aml_all_loci)

#     tbrs_hypo_bed = [] 
#     for locus in tbrs_hypo_LC.getLoci():
#         newLine = [locus.chr(), locus.start(), locus.end(), locus.ID()]
#         tbrs_hypo_bed.append(newLine)
# #    utils.unParseTable(tbrs_hypo_bed, projectFolder+'TBRS_DMR_hypo.bed', '\t')

#     tbrs_hyper_bed = []
#     for locus in tbrs_hyper_LC.getLoci():
#         newLine = [locus.chr(), locus.start(), locus.end(), locus.ID()]
#         tbrs_hyper_bed.append(newLine)
# #    utils.unParseTable(tbrs_hyper_bed, projectFolder+'TBRS_DMR_hyper.bed', '\t')



#     aml_hypo_bed = []
#     for locus in aml_hypo_LC.getLoci():
#         newLine = [locus.chr(), locus.start(), locus.end(), locus.ID()]
#         aml_hypo_bed.append(newLine)
# #    utils.unParseTable(aml_hypo_bed, projectFolder+'AML_DMR_hypo.bed', '\t')

#     aml_hyper_bed = []
#     for locus in aml_hyper_LC.getLoci():
#         newLine = [locus.chr(), locus.start(), locus.end(), locus.ID()]
#         aml_hyper_bed.append(newLine)
# #    utils.unParseTable(aml_hyper_bed, projectFolder+'AML_DMR_hyper.bed', '\t')


    # print len(mutWT_hypo_loci)
    # print len(tbrs_hypo_LC.getLoci())
    # print len(tbrs_hyper_LC.getLoci())

    # print len(aml_hypo_LC.getLoci())
    # print len(aml_hyper_LC.getLoci())


    tbrs_hypo_overlap = 0
    tbrs_hyper_overlap = 0

    aml_hypo_overlap = 0
    aml_hyper_overlap = 0

    tbrs_hypo_conserved = []

    for locus in mutWT_hypo_LC.getLoci():
#        print locus
#        overlap = tbrs_hypo_LC.getOverlap(locus, 'both')

       tbrs_hypo_overlap += len(tbrs_hypo_LC.getOverlap(locus, 'both'))
       if len(tbrs_hypo_LC.getOverlap(locus, 'both')) > 0:
           tbrs_hypo_conserved += tbrs_hypo_LC.getOverlap(locus, 'both')

       tbrs_hyper_overlap += len(tbrs_hyper_LC.getOverlap(locus, 'both'))

       aml_hypo_overlap += len(aml_hypo_LC.getOverlap(locus, 'both'))
       aml_hyper_overlap += len(aml_hyper_LC.getOverlap(locus, 'both'))

    print '%s of %s mutWT_hypo loci overlap a TBRS hypo locus' % (tbrs_hypo_overlap, len(mutWT_hypo_loci))
    print '%s of %s mutWT_hypo loci overlap a TBRS hyper locus' % (tbrs_hyper_overlap, len(mutWT_hypo_loci))
    print '==================================================================='
    print '%s of %s mutWT_hypo loci overlap a AML hypo locus' % (aml_hypo_overlap, len(mutWT_hypo_loci))
    print '%s of %s mutWT_hypo loci overlap a AML hyper locus' % (aml_hyper_overlap, len(mutWT_hypo_loci))
    print ''
    print '==================================================================='
    print len(tbrs_hypo_conserved)
    #tbrs_hypo_conserved_LC = utils.LocusCollection(tbrs_hypo_conserved)
    aml_conserved = 0
    for locus in tbrs_hypo_conserved:
        aml_overlap = aml_hypo_LC.getOverlap(locus, 'both')
        if len(aml_overlap) > 0:
            aml_conserved += len(aml_overlap)
            
    print "AML conserved", aml_conserved    

    print ''
    print '==================================================================='
    print ''

    tbrs_hypo = []    
    aml_hypo = []

    window = int(1000)
    tbrs_all_LC = utils.LocusCollection(tbrs_all_loci)
    outputFreqTable = []

    print len(mutWT_all_loci)
    print len(tbrs_all_loci)


#    mutWT_random_ctrl = random.sample(mutWT_control_loci, len(mutWT_all_loci))
#    print 'Length of randfom selection:',  len(mutWT_random_ctrl)

#     for locus in mutWT_random_ctrl:
#         #1kb
#         locusStart = locus.start()
#         locusEnd = locus.end()
#         locusChr = locus.chr()
#         locusID = locus.ID()
#         locusSense = locus.sense()

#         locusUp = locusStart - window
#         locusUpEnd = locusStart
        
#         locusDown = locusEnd
#         locusDownEnd = locusEnd + window

#         bodyBins = []
#         upstreamBins = []
#         downstreamBins = []

#         downstreamRange = list(split(range(locusDown, locusDownEnd), 50))
#         upstreamRange = list(split(range(locusUp, locusUpEnd), 50))
        
#         bodyRange = list(split(range(locusStart, locusEnd), 10))
#         for seq in upstreamRange:
#             locus = utils.Locus(locusChr, seq[0], seq[1], locusSense, '')
#             seqOverlap = tbrs_all_LC.getOverlap(locus)
#             upstreamBins.append(len(seqOverlap))

#         for seq in bodyRange:
#             try:
#                 locus = utils.Locus(locusChr, seq[0], seq[1], locusSense, '')
#                 seqOverlap = tbrs_all_LC.getOverlap(locus)
#                 bodyBins.append(len(seqOverlap))
#             except IndexError:
#                 print locus
#                 print seq

#         for seq in downstreamRange:
#             locus = utils.Locus(locusChr, seq[0], seq[1], locusSense, '')
#             seqOverlap = tbrs_all_LC.getOverlap(locus)
#             downstreamBins.append(len(seqOverlap))

        
#         outputLine = [locusID] + upstreamBins + bodyBins + downstreamBins
# #        print outputLine
#         outputFreqTable.append(outputLine)
        
#    utils.unParseTable(outputFreqTable, projectFolder+'AmishDMRvsTBRSDMR_ctrl_ratioCtrl.txt', '\t')

main()

