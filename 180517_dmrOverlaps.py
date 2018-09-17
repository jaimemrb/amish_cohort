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
#=================================FUNCTIONS======================================
#================================================================================

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


def main():

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
        if 'hypo' in line:
            tbrs_all_loci.append(utils.Locus(chrom, start, end, '.', 'tbrs_all_hypo_'+str(chrom)+':'+str(start)+'-'+str(end)))
        elif 'hyper' in line:
            tbrs_all_loci.append(utils.Locus(chrom, start, end, '.', 'tbrs_all_hyper_'+str(chrom)+':'+str(start)+'-'+str(end)))
        
    for line in aml_all:
        chrom = 'chr'+line[0]
        start = line[1]
        end = line[2]
        if 'hypo' in line:
            aml_all_loci.append(utils.Locus(chrom, start, end, '.', 'aml_all_hypo_'+str(chrom)+':'+str(start)+'-'+str(end)))
        elif 'hyper' in line:
            aml_all_loci.append(utils.Locus(chrom, start, end, '.', 'aml_all_hyper_'+str(chrom)+':'+str(start)+'-'+str(end)))

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

    print len(mutWT_hyper_loci)
    print len(mutWT_hypo_loci)

    mutWT_all_loci = mutWT_hyper_loci + mutWT_hypo_loci    
    mutWT_hypo_LC = utils.LocusCollection(mutWT_hypo_loci)

    tbrs_all_LC = utils.LocusCollection(tbrs_all_loci)
    aml_all_LC = utils.LocusCollection(aml_all_loci)
 
    tbrs_all_overlap = []
    aml_all_overlap = []

    for locus in mutWT_hypo_LC.getLoci():

        tbrs_overlap = tbrs_all_LC.getOverlap(locus, 'both')
        if len(tbrs_overlap) > 0:
            for overlapLocus in tbrs_overlap:
                overlapChrom = overlapLocus.chr()
                overlapStart = overlapLocus.start()
                overlapEnd = overlapLocus.end()

                tbrs_all_overlap.append([locus.ID(), overlapChrom, overlapStart, overlapEnd, overlapLocus.ID()])

        aml_overlap = aml_all_LC.getOverlap(locus, 'both')        
        if len(aml_overlap) > 0:
            for overlapLocus in aml_overlap:
                overlapChrom = overlapLocus.chr()
                overlapStart = overlapLocus.start()
                overlapEnd = overlapLocus.end()

                aml_all_overlap.append([locus.ID(), overlapChrom, overlapStart, overlapEnd, overlapLocus.ID()])

    utils.unParseTable(tbrs_all_overlap, projectFolder+'tables/DMRsvsTBRS_all_overlaps.txt', '\t')
    utils.unParseTable(aml_all_overlap, projectFolder+'tables/DMRsvsAML_all_overlaps.txt', '\t')

main()


