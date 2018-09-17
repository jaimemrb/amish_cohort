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
    
    genome = 'hg19'

    projectFolder = '/storage/goodell/projects/jmreyes/amish_ayala/'
    bedFolder = projectFolder+'bed/'
    canyonBed = bedFolder+'canyon_Mut_sizeSelected.bed'
    scripts = projectFolder+'scripts/'
    sampleName = 'Mut'
    outputUp = 'HG19_'+sampleName+'_1000extend_upstreamFlanking.bed'
    outputDown = 'HG19_'+sampleName+'_'+'1000extend_downstreamFlanking.bed'

    geneGTF = '/storage/goodell/home/jmreyes/grail/genomes/Homo_sapiens/UCSC/%s/Annotation/Genes/genes.gtf'%(genome)

    bedList = [outputUp, outputDown]
    cmdBash = [['#!/usr/bin/bash']]
    cmdOut = projectFolder+'scripts/'+sampleName+'_geneIntersect.sh'
    for bed in bedList:
        bedIn = bed
        bedName = bedIn.split('/')[-1].split('.bed')[0]
        sortedOut = bedFolder+bedName+'.sorted.bed'
        intersectOut = bedFolder+bedName+'_geneIntersect.bed'
        sortBedCmd = 'sort -k1,1 -k2,2n %s > %s' % (bedFolder+bedIn, sortedOut)
        cmdBash.append([sortBedCmd])
        
        intersectCmd = 'bedtools closest -d -a %s -b %s > %s' % (sortedOut, geneGTF, intersectOut)
        cmdBash.append([intersectCmd])

    utils.unParseTable(cmdBash, cmdOut, '\t')

main()

