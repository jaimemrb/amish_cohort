#!/usr/bin/python

#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================

import sys

print "Using python version %s" % sys.version

sys.path.append('/storage/goodell/home/jmreyes/pipeline_jr/')
import utils
import pipeline_dfci


#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================


#================================================================================
#===================================CLASSES======================================
#================================================================================


#================================================================================
#=================================FUNCTIONS======================================
#================================================================================

def expansionStat(locus1, locus2, expansion = 0.1):
    '''uses locus as list [chr, start, end, length]'''
    locusLength1 = int(locus1[3])
    locusLength2 = int(locus2[3])
    expandedList = []

    if locusLength1 - locusLength2 > expansion * locusLength1:        
#        print locus1, locus2, locusLength1, locusLength2
        newLine = locus2
    else:
        newLine = []
    return newLine

def getExpanded(locusTable, expansion, status, output):
    
    loci = utils.parseTable(locusTable, '\t')
    expandedList = []
    for line in loci:
        wtLocus = line[0:4]
        mutLocus = line[4:8]
        if status == 'WT':
            newLine = expansionStat(wtLocus, mutLocus, expansion = 0.1)
            if len(newLine) >0:
                expandedList.append(newLine)
        elif status == 'MUT':
            newLine = expansionStat(mutLocus, wtLocus, expansion = 0.1)
            if len(newLine) > 0:
                expandedList.append(newLine)

    print len(expandedList), ' expanded loci in ', status
    utils.unParseTable(expandedList, output, '\t')

        
#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


def main():

    projectFolder = '/storage/goodell/projects/jmreyes/amish_ayala/'
    bedFolder = projectFolder+'bed/'
    wtExpanded = projectFolder+'tables/WT_canyons_expanded.txt'
    mutExpanded = projectFolder+'tables/MUT_canyons_expanded.txt'

    getExpanded(mutExpanded, expansion = 0.1, status = 'MUT', output = projectFolder+'bed/MUT_canyon_expansion_01.bed')
    getExpanded(wtExpanded, expansion = 0.1, status = 'WT', output = projectFolder+'bed/WT_canyon_expansion_01.bed')
    




main()

