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




#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


def main():

    projectFolder = '/storage/goodell/projects/jmreyes/amish_ayala/'
    bedFolder = projectFolder+'bed/'
    wtCanyonBed = bedFolder+'canyon_WT_sizeSelected.bed'
    mutCanyonBed = bedFolder+'canyon_Mut_sizeSelected.bed'
    
    wtCanyonLocusCollection = utils.LocusCollection([utils.Locus(x[0], x[1], x[2], '.', 'wt_'+str(x[0])+':'+str(x[1])+'-'+str(x[2])) for x in utils.parseTable(wtCanyonBed, '\t')])
    mutCanyonLocusCollection = utils.LocusCollection([utils.Locus(x[0], x[1], x[2], '.', 'mut_'+str(x[0])+':'+str(x[1])+'-'+str(x[2])) for x in utils.parseTable(mutCanyonBed, '\t')])
    overlappingCanyons = []

    wtExpansion = []
    mutExpansion = []

    wtUnique = []
    mutUnique = []
    overlapCounter = 0
    mutOverlap = 0
    for locus in wtCanyonLocusCollection.getLoci():
        wtMutOverlap = mutCanyonLocusCollection.getOverlap(locus, 'both')
        if len(wtMutOverlap) > 0:
            overlapCounter += 1
            for overlap in wtMutOverlap:
                newLine = [locus.chr(), locus.start(), locus.end(), locus.end()-locus.start(), overlap.chr(), overlap.start(), overlap.end(), overlap.end()-overlap.start()]
                wtLength = locus.end()-locus.start()
                mutLength = overlap.end()-overlap.start()
                if mutLength > wtLength:
                    mutExpansion.append(newLine)
                elif wtLength > mutLength:
                    wtExpansion.append(newLine)
        else:
            wtUnique.append(locus)

    for locus in mutCanyonLocusCollection.getLoci():
        mutWTOverlap = wtCanyonLocusCollection.getOverlap(locus, 'both')
        if len(mutWTOverlap) > 0:
            mutOverlap += 1
        else:
            mutUnique.append(locus)


    print len(mutExpansion)
    print len(wtExpansion)
    utils.unParseTable(mutExpansion, projectFolder+'tables/MUT_canyons_expanded.txt', '\t')
    utils.unParseTable(wtExpansion, projectFolder+'tables/WT_canyons_expanded.txt', '\t')

main()

