
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



#================================================================================
#============================GLOBAL PARAMETERS===================================
#================================================================================

#add locations of files and global parameters in this section


dataFile ='/location/file.txt'
genome = 'hg19'
projectFolder = '/storage/goodell/projects/jmreyes/amish_ayala/'

#================================================================================
#===================================CLASSES======================================
#================================================================================

#user defined classes here

#================================================================================
#=================================FUNCTIONS======================================
#================================================================================

#write your specific functions here
annotFile = '/storage/goodell/home/jmreyes/pipeline/annotation/%s_refseq.ucsc' % (genome)

startDict = utils.makeStartDict(annotFile)
startLoci = []
#for TR, -30, +300 and genebody +0

for gene in startDict.keys():
    geneChrom = startDict[gene]['chr']
    geneStart = startDict[gene]['start']
    geneEnd   = startDict[gene]['end']
    geneSense  = startDict[gene]['sense']
    
#    newLocus  = [geneChrom, gene, '', geneStart]

    
    newLocus = utils.makeTSSLocus(gene, startDict, 0, 0)
    

    startLoci.append([newLocus.chr(), newLocus.start(), newLocus.end()+1, newLocus.sense(), startDict[gene]['name'], newLocus.ID()])

utils.unParseTable(startLoci, projectFolder+'HG19_genes.bed', '\t')




 
#================================================================================
#===============================MAIN RUN=========================================
#================================================================================


#write the actual script here




# def main():

#  #source gene list
#  leyTable = utils.parseTable(projectFolder+'total.rnaseq.counts.tsv', '\t')
 
#  ensemblENST = utils.parseTable(projectFolder+'Homo_sapiens.GRCh37.85.ena.tsv', '\t')
#  enstToGene = utils.parseTable(projectFolder+'hg19_ensemblToGene', '\t')
#  ensgToEnstDict = {}
 
#  for line in ensemblENST[1:]:
#   ensg = line[2]
#   enst = line[3]
#   ensgToEnstDict[ensg] = enst
  
#  enstToGeneDict = {}
#  for line in enstToGene[1:]:
#   enst = line[0]
#   gene = line[1]
#   enstToGeneDict[enst] = gene
  
#  header = ['GENE', 'ENST_ID', 'ENSG_ID'] + leyTable[0]
 
#  newTable = [header]

#  for line in leyTable[1:]:
#   ensemblID = line[0].split('.')[0]
#   if ensemblID in ensgToEnstDict.keys():
#    enstName = ensgToEnstDict[ensemblID]
#    geneName = enstToGeneDict[enstName]
#   else:
#    enstName = 'NA'
#    geneName = 'NA'

#   newLine = [geneName, enstName] + line
#   newTable.append(newLine)

#  utils.unParseTable(newTable, projectFolder+'Ley_geneNames_rnaCounts.tsv', '\t')
# main()


