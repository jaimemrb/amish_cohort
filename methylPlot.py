 #!/usr/bin/python


#================================================================================
#=============================DEPENDENCIES=======================================
#================================================================================

import sys

print "Using python version %s" % sys.version


#importing utils package
sys.path.append('/storage/goodell/home/jmreyes/pipeline_jr/')
import utils



#================================================================================
#=================================FUNCTIONS======================================
#================================================================================


#split regions into UP, BODY, DOWN
def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

def writeSplitBeds(bed, analysisName, outputFolder, window=50):

    dmrList = [x for x in utils.parseTable(bed, '\t')]
#    print len(refGenes)

    dmrBed = []
    endsBed = []
    startsBed = []

    for line in dmrList:
        dmrID = line[3]

        dmrCoords = [line[0], int(line[1]), int(line[2]), dmrID]

        dmrBed.append(dmrCoords)

        startExtend =  [line[0], int(line[1])-window, int(line[1]), dmrID]
        endExtend =  [line[0], int(line[2]), int(line[2])+window,dmrID]
        
        endsBed.append(endExtend)
        startsBed.append(startExtend)

    print len(dmrBed)
    utils.unParseTable(exonBed, outputFolder+'bed/'+analysisName+'_BODY_-0_+0.bed', '\t')
    print len(startsBed)
    utils.unParseTable(startsBed, outputFolder+'bed/'+analysisName+'_UPSTREAM_-'+str(window)+'_+'+str(window)+'.bed', '\t')
    print len(endsBed)
    utils.unParseTable(endsBed, outputFolder+'bed/'+analysisName+'_DOWNSTREAM_-'+str(window)+'_+'+str(window)+'.bed', '\t')

    return [outputFolder+'bed/'+analysisName+'_BODY_-0_+0.bed', outputFolder+'bed/'+analysisName+'_UPSTREAM_-'+str(window)+'_+0.bed', outputFolder+'bed/'+analysisName+'_DOWNSTREAM_-0_+'+str(window)+'.bed']

#getCGOverlaps
#def getCpgOverlap(bed, wgbs):



#================================================================================
#===============================MAIN RUN=========================================
#================================================================================

def main():

    from optparse import OptionParser
    usage = "usage: %prog [options] -i [INPUTFILE] -o [OUTPUT]"
    parser = OptionParser(usage=usage)    
    parser.add_option("-i", "--i", dest="input", nargs=1, default=None,
                      help="Enter path to WGBS calls")
    parser.add_option("-b", "--b", dest="bed", nargs=1, default=None,
                      help="Enter path to DMR bed")
    parser.add_option("-n", "--n", dest="binNum", nargs=1, default=50,
                      help="Enter bin number")
    parser.add_option("-o", "--o", dest="output", nargs=1, default=None,
                      help="Enter path to output table, will be tab-delim by default")

    (options, args) = parser.parse_args()

    if not options.input or not options.output:
        parser.print_help()
        exit()

#    projectFolder = '/storage/goodell/home/jmreyes/projects/amish_ayala/'

    #import data table
    print('Loading table')

    
    #import regions
    dmrTable = utils.parseTable(options.bed, '\t')
    dmrBed = sorted(dmrTable, key=lambda x:(x[0], int(x[1]), int(x[2])))
    print dmrBed[1:10]
    chromList = utils.uniquify([x[0] for x in dmrBed if 'chr' in x[0]])
    print chromList

    binNum = int(options.binNum)

    print('Making Dict')
    methylCallsDict = {}
    for chrom in chromList:
        #load chromosome CpGs
        methylCallsDict[chrom] = {}
        with open(options.input) as f:
            results = [line for line in f if line.startswith(chrom)]
        ticker = 0
        for x in results:
            ticker += 1
            if ticker % 1000000 == 0:
                print ticker, ' CpGs loaded to Dict'
            line = x.split('\t')
            cpgStart = line[1]
            ratio = line[4]
            methylCallsDict[chrom][cpgStart] = ratio
        f.close()


    print("Checking CpGs in regions")
    methylTable = open(options.output, 'w')
    ticker = 0
    for line in dmrBed:
        ticker += 1
        if ticker % 100 == 0:
            print ticker, " DMRs mapped"
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        locusID = line[4]
        if end-start > binNum:
            dmrSeq = list(split(range(start, end), binNum))
        else:
            shortRange = range(start, end)
            diff = binNum - len(shortRange)
            newEnd = end + diff
            dmrSeq = list(split(range(start, newEnd), binNum))
            

        methylLine = []
        keyList = [x for x in methylCallsDict[chrom].keys() if start <= int(x) <= end]
        
        cpgCount = 0
        
        for i in dmrSeq:
            try:
                binMin = min(i)
                binMax = max(i)
            except ValueError:
                print i
                print dmrSeq
            dmrCpgs = [float(methylCallsDict[chrom][x]) for x in keyList if binMin <= int(x) <= binMax]
            if len(dmrCpgs) > 0:
                dmrMean = sum(dmrCpgs)/len(dmrCpgs)
                cpgCount += len(dmrCpgs)
            else:
                dmrMean = 'NA'
                cpgCount += len(dmrCpgs)

            methylLine.append(dmrMean)
        outLine = [locusID, str(cpgCount)] + methylLine
        methylTable.write('\t'.join([str(x) for x in outLine]) + '\n') 

    print("Complete")
    methylTable.close()
#    utils.unParseTable(methylTable, options.output, '\t')


if __name__ == "__main__":
    main()

