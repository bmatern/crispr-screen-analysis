import argparse
from os import listdir, makedirs
from os.path import isfile, join, isdir

import gzip
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

import multiprocessing as mp


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--sequencedlibrary", help="directory containing fastq or fastq.gz files", required=True)
    parser.add_argument("-f", "--fasta", help="directory to write blast input fasta files", required=True)
    parser.add_argument("-br", "--blastresults", help="directory to write blast input fasta files", required=True)
    parser.add_argument("-bd", "--blastdatabase", help="blast database identifier", required=True)
    parser.add_argument("-s", "--summarydirectory", help="summary results directory", required=True)
    parser.add_argument("-t", "--threads", help="threadcount for blast analysis", required=False, default=1)
    parser.add_argument("-bs", "--batch", help="batch size", required=False, default=1000)
    args = parser.parse_args()

    return args


def loadReadFastqFiles(libraryDirectory=None):
    print('Loading read files from directory:' + str(libraryDirectory))

    readData = {}
    readCount = 0

    fileList = [f for f in listdir(libraryDirectory) if isfile(join(libraryDirectory, f))]

    for fileName in fileList:
        readData[fileName] = {}
        # If it's a fastq.gz file
        if(fileName.lower().endswith('.fastq.gz')):
            with gzip.open(join(libraryDirectory, fileName), "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    readCount += 1
                    currentRead = {}
                    currentRead['seq'] = str(record.seq)
                    currentRead['qual'] = record.letter_annotations['phred_quality']
                    readData[fileName][str(record.id)] = currentRead

        else:
            raise Exception ('Have not implemented handling this file type:' + str(fileName))

    print ('Read ' + str(readCount) + ' reads from ' + str(len(fileList)) + ' files.')
    return readData


def writeReadFastas(readData=None, outputDirectory=None, batchSize=None, newline='\n'):
    print('Writing Read Fastas to directory:' + str(outputDirectory))

    if int(batchSize) < 1:
        raise Exception('Invalid batch size:' + str(batchSize))

    if not isdir(outputDirectory):
        makedirs(outputDirectory)


    currentReadCount=0
    currentBatchIndex=0
    fileHandle=None

    for readFileName in readData.keys():
        for readId in readData[readFileName].keys():
            currentReadCount += 1

            if fileHandle is None:
                currentBatchIndex += 1
                outputFileName = join(outputDirectory, 'Batch_' + str("{:04d}".format(currentBatchIndex)) + '.fasta')
                fileHandle = open(outputFileName, 'w')

            fileHandle.write('>' + str(readId) + newline)
            fileHandle.write(str(readData[readFileName][readId]['seq']) + newline)

            if currentReadCount >= batchSize:
                fileHandle.close()
                fileHandle = None
                currentReadCount = 0


    if fileHandle is not None:
        fileHandle.close()


def parseBlastXml(xmlFileName=None, resultFileName=None, deleteXmlFile=False, delimiter=',', newline='\n'):
    with open(resultFileName, 'w') as resultFile:
        blastRecords = NCBIXML.parse(open(xmlFileName, 'r'))

        for blastRecord in blastRecords:
            numberAlignments = len(blastRecord.alignments)
            queryName = str(blastRecord.query)

            if numberAlignments<=1 :
                resultFile.write(queryName + delimiter + newline)
            else:
                bestHit = blastRecord.alignments[0].title.replace(delimiter,';')
                resultFile.write(queryName + delimiter + bestHit + newline)

    # TODO: Delete the XML File.
    # Delete XML file.


def performBlastBatch(fastaDirectory=None, fastaFile=None, outputDirectory=None, blastDatabase=None):
    print('Blasting file ' + str(fastaFile))
    outputFile = join(outputDirectory,fastaFile + ".xml")
    cline = NcbiblastnCommandline(query=join(fastaDirectory,fastaFile), db=blastDatabase, out=outputFile, outfmt=5)
    cline()

    parseBlastXml(xmlFileName=outputFile, resultFileName=outputFile.replace('.fasta.xml','.csv'), deleteXmlFile=True)
    print('Done Blasting file ' + str(fastaFile))



def blastReads(fastaDirectory=None, outputDirectory=None, threadCount=1, blastDatabase=None):
    fastaFileList = sorted([f for f in listdir(fastaDirectory) if (isfile(join(fastaDirectory, f)) and f.endswith('.fasta'))])
    print('Running Blast on ' + str(len(fastaFileList)) + ' fasta files, using ' + str(threadCount) + ' threads')

    if not isdir(outputDirectory):
        makedirs(outputDirectory)

    pool = mp.Pool(threadCount)

    for fastaFile in fastaFileList:
        if(threadCount <= 1):
            performBlastBatch(fastaDirectory=fastaDirectory, fastaFile=fastaFile, outputDirectory=outputDirectory, blastDatabase=blastDatabase)
        else:
            pool.starmap_async(performBlastBatch, [[fastaDirectory, fastaFile, outputDirectory, blastDatabase]])

    pool.close()
    pool.join()


def findOriginReadFile(readData=None, readName=None):
    for fileName in readData:
        if readName in readData[fileName].keys():
            return fileName
    return '?'


def summarizeAndCombineBlastResults(readData=None, blastResultDirectory=None, finalResultDirectory=None, newline='\n', delimiter=','):
    resultFileList = sorted([f for f in listdir(blastResultDirectory) if (isfile(join(blastResultDirectory, f)) and f.endswith('.csv'))])
    print('Combining ' + str(len(resultFileList)) + ' blast results in :' + str(blastResultDirectory))

    perReadFileBlastHits = {}
    totalBlastHits = {}

    outputCombinedFileName = join(finalResultDirectory,'CombinedBlastHits.csv')
    with open(outputCombinedFileName,'w') as outputCombinedFile:
        for resultFileName in resultFileList:
            with open(join(blastResultDirectory, resultFileName), 'r') as resultFile:
                resultFileLines = resultFile.readlines()
                for resultFileLine in resultFileLines:
                    readName, bestHit = resultFileLine.strip().split(delimiter)
                    originReadFile = findOriginReadFile(readData=readData, readName=readName)

                    outputCombinedFile.write(delimiter.join([originReadFile,readName,bestHit]) + newline)

                    # Summary data
                    if bestHit is not None and len(bestHit.strip())>1:
                        # Total summary
                        if bestHit not in totalBlastHits.keys():
                            totalBlastHits[bestHit]=0
                        totalBlastHits[bestHit] += 1

                        # Per-Readfile summary
                        if originReadFile not in perReadFileBlastHits.keys():
                            perReadFileBlastHits[originReadFile] = {}
                        if bestHit not in perReadFileBlastHits[originReadFile].keys():
                            perReadFileBlastHits[originReadFile][bestHit] = 0
                        perReadFileBlastHits[originReadFile][bestHit] += 1


    summaryFileName = join(finalResultDirectory,'HitSummary.txt')
    with open(summaryFileName,'w') as summaryFile:
        summaryFile.write('All reads blast hits:' + newline + newline)
        sortedHits = {k: v for k, v in sorted(totalBlastHits.items(), key=lambda item: item[1], reverse=True)}
        for bestHit in sortedHits.keys():
            summaryFile.write(bestHit + '  (' + str(sortedHits[bestHit]) + ' reads)' + newline)

        summaryFile.write(newline + newline)

        for readFileName in sorted(list(perReadFileBlastHits.keys())):
            summaryFile.write(readFileName + ' blast hits:' + newline + newline)

            sortedHits = {k: v for k, v in sorted(perReadFileBlastHits[readFileName].items(), key=lambda item: item[1], reverse=True)}
            for bestHit in sortedHits.keys():
                summaryFile.write(bestHit + '  (' + str(sortedHits[bestHit]) + ' reads)' + newline)

            summaryFile.write(newline + newline)




if __name__ == '__main__':
    args = parseArgs()

    # TODO: Maybe need to combine this, and write to hard drive instead of keeping it all in memory.
    readData = loadReadFastqFiles(libraryDirectory=args.sequencedlibrary)
    writeReadFastas(readData=readData, outputDirectory=args.fasta, batchSize=int(args.batch))
    blastReads(fastaDirectory=args.fasta, outputDirectory=args.blastresults, blastDatabase=args.blastdatabase, threadCount=int(args.threads))
    summarizeAndCombineBlastResults(readData=readData, blastResultDirectory=args.blastresults, finalResultDirectory=args.summarydirectory)

    print('Done.')