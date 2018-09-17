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
    canyonBed = bedFolder+'canyon_WT_sizeSelected.bed'

    extension = 1000
    sampleName = 'WT'
    upstreamEdgeBed = []
    downstreamEdgeBed = []

    outputUp = 'HG19_'+sampleName+'_'+str(extension)+'extend_upstreamFlanking.bed'
    outputDown = 'HG19_'+sampleName+'_'+str(extension)+'extend_downstreamFlanking.bed'
    for line in utils.parseTable(canyonBed, '\t'):
        print line
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        canyon_name = sampleName+'_'+str(chrom)+'(.):'+str(start)+'-'+str(end)
        
        if start > extension:
            startUpstream = start - extension
            startDownstream = start + extension

            endUpstream = end - extension
            endDownstream = end + extension
            upstreamLine = [chrom, startUpstream, startDownstream, canyon_name+'_5Flank']
            downstreamLine = [chrom, endUpstream, endDownstream, canyon_name+'_3Flank']
            upstreamEdgeBed.append(upstreamLine)
            downstreamEdgeBed.append(downstreamLine)

        else:
            pass



    utils.unParseTable(upstreamEdgeBed, bedFolder+outputUp, '\t')
    utils.unParseTable(downstreamEdgeBed, bedFolder+outputDown, '\t')
        

main()

