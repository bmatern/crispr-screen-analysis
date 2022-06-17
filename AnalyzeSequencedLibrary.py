import argparse
from os import listdir, makedirs
from os.path import isfile, join, isdir

import gzip
from traceback import print_exc

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="darkgrid")

import multiprocessing as mp

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--sequencedlibrary", help="directory containing fastq or fastq.gz files", required=True)
    parser.add_argument("-f", "--fasta", help="directory to write blast input fasta files", required=True)
    parser.add_argument("-br", "--blastresults", help="directory to write blast input fasta files", required=True)
    parser.add_argument("-bd", "--blastdatabase", help="blast database identifier", required=True)
    parser.add_argument("-s", "--summarydirectory", help="summary results directory", required=True)
    parser.add_argument("-g", "--guidernas", help="guide RNA Sequence file", required=True)
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
                    #if bestHit is not None and len(bestHit.strip())>1:
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
            if bestHit is None or len(bestHit.strip())<1:
                hitName='No Hits Found '
            else:
                hitName=bestHit
            summaryFile.write(hitName + '  (' + str(sortedHits[bestHit]) + ' reads)' + newline)

        summaryFile.write(newline + newline)

        for readFileName in sorted(list(perReadFileBlastHits.keys())):
            summaryFile.write(readFileName + ' blast hits:' + newline + newline)

            sortedHits = {k: v for k, v in sorted(perReadFileBlastHits[readFileName].items(), key=lambda item: item[1], reverse=True)}
            for bestHit in sortedHits.keys():
                if bestHit is None or len(bestHit.strip()) < 1:
                    hitName = 'No Hits Found '
                else:
                    hitName = bestHit
                summaryFile.write(hitName + '  (' + str(sortedHits[bestHit]) + ' reads)' + newline)

            summaryFile.write(newline + newline)

def writeReadSummary(readData=None, finalResultDirectory=None, newline='\n',delimiter=','):
    print('Writing Read Summary to ' + str(finalResultDirectory))

    if not isdir(finalResultDirectory):
        makedirs(finalResultDirectory)

    readLengths=[]
    readLengthsPerFile={}

    readLengthFilename = join(finalResultDirectory,'ReadLengths.csv')
    with open (readLengthFilename,'w') as readLengthFile:
        for readFileName in readData.keys():
            readLengthsPerFile[readFileName]=[]
            for readName in readData[readFileName].keys():
                # TODO: Calculate average read quality, maybe
                readLengthFile.write(delimiter.join([readFileName,readName,str(len(readData[readFileName][readName]['seq']))]) + newline)

                readLengths.append(len(readData[readFileName][readName]['seq']))
                readLengthsPerFile[readFileName].append(len(readData[readFileName][readName]['seq']))


    sns.kdeplot(readLengths)
    plt.title('Read Length, All ' + ' (N=' + str(len(readLengths)) + ')')
    plt.savefig(join(finalResultDirectory,'Length.AllReads.png'), bbox_inches='tight')
    plt.clf()
    for readFileName in readLengthsPerFile.keys():
        sns.kdeplot(readLengthsPerFile[readFileName])
        plt.title('Read Length, ' + str(readFileName) + ' (N=' + str(len(readLengthsPerFile[readFileName])) + ')')
        plt.savefig(join(finalResultDirectory,'Length.' + str(readFileName) + '.png'), bbox_inches='tight')
        plt.clf()

def readGuideRNAs(guideFileName=None, delimiter='\t', header=True):
    print('Reading Guide RNA File:' + str(guideFileName))

    guideRNAs = {}
    with open(guideFileName,'r') as guideFile:
        dataLines = guideFile.readlines()
        for lineIndex, line in enumerate(dataLines):
            if not header or lineIndex > 0:
                dataTokens = line.split(delimiter)

                geneName = dataTokens[1]
                sgRnaTarget = dataTokens[6]

                if geneName not in guideRNAs.keys():
                    guideRNAs[geneName] = []
                guideRNAs[geneName].append(sgRnaTarget)

    return guideRNAs

def checkReadsAgainstGuideRNAs(readData=None, guideRnas=None, threadCount=1, finalResultDirectory=None, delimiter=',', newline='\n'):
    print('Checking Reads against Guide RNAs, threadcount=' + str(threadCount))

    pool = mp.Pool(threadCount)
    results = []

    for readfileIndex, readFileName in enumerate(sorted(list(readData.keys()))):
        print('Checking Read File ' + str(readFileName) + ' (' + str(readfileIndex+1) + ' of ' + str(len(list(readData.keys()))) + ')')

        #compareReadFileAgainstGuideRNAs(readFileName=readFileName, reads=readData[readFileName], guideRnas=guideRnas)
        results.append(pool.starmap_async(compareReadFileAgainstGuideRNAs, [[readFileName, readData[readFileName], guideRnas]]))

    print('Done checking Reads against Guide RNAs')

    pool.close()
    pool.join()

    combinedHitSummary = {}
    combinedReadSummary = {}

    for resultObject in results:
        hitResult, readResult = resultObject.get()[0]

        #print('I found a hit result with ' + str(len(hitResult.keys())) + ' hits.')
        #print('I found a read result with ' + str(len(readResult.keys())) + ' reads.')

        for geneName in hitResult.keys():
            if geneName not in combinedHitSummary.keys():
                combinedHitSummary[geneName] = set()
            combinedHitSummary[geneName] = combinedHitSummary[geneName].union(hitResult[geneName])

        for readName in readResult.keys():
            combinedReadSummary[readName] = readResult[readName]

    # Sort by read counts
    combinedHitSummary = {k: v for k, v in sorted(combinedHitSummary.items(), key=lambda item: len(list(item[1])), reverse=True)}

    hitSummaryFileName = join(finalResultDirectory, 'sgRNAHitSummary.csv')
    with open(hitSummaryFileName,'w') as hitSummaryFile:
        hitSummaryFile.write(delimiter.join(['Total Genes Found:' , str(len(list(combinedHitSummary.keys())))]) + newline)

        hitSummaryFile.write(newline+newline + delimiter.join(['Gene','Read Count']) + newline)


        for geneName in combinedHitSummary.keys():
            hitSummaryFile.write(delimiter.join([geneName, str(len(list(combinedHitSummary[geneName])))]) + newline)

        for geneName in combinedHitSummary.keys():
            hitSummaryFile.write(newline + newline + delimiter.join([geneName, str(len(list(combinedHitSummary[geneName])))]) + newline)
            for readName in list(combinedHitSummary[geneName]):
                hitSummaryFile.write(delimiter.join([findOriginReadFile(readData=readData, readName=readName),readName]) + newline)

    perReadSummaryFileName = join(finalResultDirectory, 'readSummary.csv')
    with open(perReadSummaryFileName,'w') as readSummaryFile:
        readSummaryFile.write(delimiter.join(['Read File', 'Read Name', 'Read Length', 'Guide RNA Gene', 'Guide RNA Sequence']) + newline)
        for readName in combinedReadSummary.keys():
            readFileName = findOriginReadFile(readData=readData, readName=readName)
            readSummaryFile.write(delimiter.join([readFileName, readName
                , str(combinedReadSummary[readName]['read_length'])
                , combinedReadSummary[readName]['sg_gene_found']
                , combinedReadSummary[readName]['sg_rna_found']]) + newline)


def replaceText(fullSeq=None, replaceText=None, startIndex=None):
    temp = list(fullSeq)
    for i in range(0, len(replaceText)):
        temp[i+startIndex] = replaceText[i]
    returnText = "".join(temp)
    return returnText


def writeReadAlignmentToConsole(readName=None, readSummary=None, readSeq=None, newline= '\n'):
    print('Hit! Readname:' + str(readName) + ', Genename:' + str(readSummary['sg_gene_found']))
    print(readSeq)

    alignmentSeq = ' ' * len(readSeq)
    # Aligned sgRNA Sequence
    alignmentSeq = replaceText(fullSeq=alignmentSeq, replaceText=readSummary['sg_rna_found'], startIndex = readSeq.index(readSummary['sg_rna_found']))

    # Aligned Flanking Sequence
    alignmentSeq = replaceText(fullSeq=alignmentSeq, replaceText=readSummary['gecko_flank_left_seq'].lower(), startIndex=readSummary['gecko_flank_left_index'])
    alignmentSeq = replaceText(fullSeq=alignmentSeq, replaceText=readSummary['gecko_flank_right_seq'].lower(), startIndex=readSummary['gecko_flank_right_index'])

    alignmentSeq = replaceText(fullSeq=alignmentSeq, replaceText=readSummary['gecko_reverse_primer_seq'], startIndex=readSummary['gecko_reverse_primer_index'])

    #['gecko_reverse_primer_seq'], readSummary[readName]['gecko_reverse_primer_index']

    print(alignmentSeq + newline)

    pass

def findSubsequenceInRead(readSeq=None, subSequences=None, subsequenceMinimum=12):
    for subSequence in subSequences:
        # TODO: This will break if subsequence < 10
        for trimIndex in range(0,len(subSequence)-subsequenceMinimum):
            # trim from left
            subsequenceTrimmed = subSequence[trimIndex:]
            if subsequenceTrimmed in readSeq:
                return subsequenceTrimmed, readSeq.index(subsequenceTrimmed)

            # trim from right
            subsequenceTrimmed = subSequence[0:len(subSequence)-trimIndex]
            if subsequenceTrimmed in readSeq:
                return subsequenceTrimmed, readSeq.index(subsequenceTrimmed)

    return '',0


def compareReadFileAgainstGuideRNAs(readFileName=None, reads=None, guideRnas=None):
    print('Comparing reads for file ' + str(readFileName))

    # Reverse Complement Guide RNA Sequences
    for geneName in guideRnas.keys():
        guideRnas[geneName].extend([str(Seq(sgRna).reverse_complement()) for sgRna in guideRnas[geneName]])

    # Gecko Flanking Sequences
    # https://github.com/fengzhanglab/Screening_Protocols_manuscript/blob/master/design_library.py
    geckoFlankLeft = ['TTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG']
    geckoFlankLeft.extend([str(Seq(geckoSeq).reverse_complement()) for geckoSeq in geckoFlankLeft])
    geckoFlankRight = ['GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGT']
    geckoFlankRight.extend([str(Seq(geckoSeq).reverse_complement()) for geckoSeq in geckoFlankRight])

    # Some reverse primer sequences, I think i see these in the reads
    # https://www.nature.com/articles/nprot.2017.016/tables/3
    geckoLibraryReverseSequence = ['GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCGACTCGGTGCCACTTTTTCAA']
    geckoLibraryReverseSequence.extend([str(Seq(geckoSeq).reverse_complement()) for geckoSeq in geckoLibraryReverseSequence])

    guideRnaHitSummary = {}
    # TODO: I'm using read name as a key in this dictionary. That's flawed because the same read name can appear in multiple files.
    readSummary={}

    for readName in reads.keys():

        try:


            readSeq = reads[readName]['seq']

            readSummary[readName] = {}

            #readSummary[readName]['direction'] = None
            readSummary[readName]['sg_rna_found']=''
            readSummary[readName]['sg_gene_found']=''
            readSummary[readName]['read_length']=len(readSeq)

            # Check "forward" target sequences
            for geneName in guideRnas.keys():
                for sgRna in guideRnas[geneName]:
                    if sgRna in readSeq:
                        readSummary[readName]['sg_rna_found']=sgRna
                        readSummary[readName]['sg_gene_found']=geneName
                        #readSummary[readName]['direction']='Forward'

                        if geneName not in guideRnaHitSummary.keys():
                            guideRnaHitSummary[geneName] = set()

                        guideRnaHitSummary[geneName].add(readName)

                        # We're only expecting one guide RNA per read....right?
                        break

            # Other sequences found in the reads...
            readSummary[readName]['gecko_flank_left_seq'], readSummary[readName]['gecko_flank_left_index'] = findSubsequenceInRead(readSeq=readSeq, subSequences=geckoFlankLeft)
            readSummary[readName]['gecko_flank_right_seq'], readSummary[readName]['gecko_flank_right_index'] = findSubsequenceInRead(readSeq=readSeq, subSequences=geckoFlankRight)

            readSummary[readName]['gecko_reverse_primer_seq'], readSummary[readName]['gecko_reverse_primer_index'] = findSubsequenceInRead(readSeq=readSeq, subSequences=geckoLibraryReverseSequence)

            if(readSummary[readName]['sg_rna_found'] is not None):
                writeReadAlignmentToConsole(readName=readName, readSummary=readSummary[readName],  readSeq=readSeq)

        except Exception as e:
            print('EXCEPTION in read ' + str(readName) + str(e))
            print_exc()
            raise(e)

    print('Done comparing reads for file ' + str(readFileName))
    return guideRnaHitSummary, readSummary

if __name__ == '__main__':
    args = parseArgs()

    # TODO: Maybe need to combine this, and write to hard drive instead of keeping it all in memory.
    readData = loadReadFastqFiles(libraryDirectory=args.sequencedlibrary)
    #writeReadFastas(readData=readData, outputDirectory=args.fasta, batchSize=int(args.batch))
    #writeReadSummary(readData=readData, finalResultDirectory=args.summarydirectory)
    guideRnas = readGuideRNAs(guideFileName=args.guidernas)
    checkReadsAgainstGuideRNAs(readData=readData, guideRnas=guideRnas, threadCount=int(args.threads), finalResultDirectory=args.summarydirectory)
    #blastReads(fastaDirectory=args.fasta, outputDirectory=args.blastresults, blastDatabase=args.blastdatabase, threadCount=int(args.threads))
    #summarizeAndCombineBlastResults(readData=readData, blastResultDirectory=args.blastresults, finalResultDirectory=args.summarydirectory)

    print('Done.')