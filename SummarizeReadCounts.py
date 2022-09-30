import argparse
from os import listdir, makedirs
from os.path import isfile, join, isdir


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--results", help="directory with result files from count_spacers.py", required=True)
    parser.add_argument("-g", "--guidernas", help="guide RNA Sequence file", required=True)
    args = parser.parse_args()

    return args


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


def analyzeDirectory(resultDirectory=None, guideRnas=None, delimiter=',', newline='\n'):
    print('Searching for result files in ' + str(resultDirectory))

    libraryCountFileNames = sorted([f for f in listdir(resultDirectory) if (isfile(join(resultDirectory, f)) and f.endswith('.library_count.csv'))])
    statisticsFileNames = sorted([f for f in listdir(resultDirectory) if (isfile(join(resultDirectory, f)) and f.endswith('.statistics.txt'))])

    # Read Statistics Files
    statisticNames = set()
    statisticsData = {}
    for statisticsFileName in statisticsFileNames:
        print('Reading statistics file ' + str(statisticsFileName))
        key = statisticsFileName.replace('.statistics.txt','')
        statisticsData[key] = {}
        with open(join(resultDirectory, statisticsFileName),'r') as statisticsFile:
            for line in list(statisticsFile.readlines()):
                statisticName, value = line.strip().split(':')
                statisticsData[key][statisticName.strip()] = value.strip()
                statisticNames.add(statisticName)

    # Write Statistics Summary
    statsSummaryFileName = join(resultDirectory, 'StatsSummary.csv')
    with open(statsSummaryFileName,'w') as statsSummaryFile:
        datasetNames = list(statisticsData.keys())
        statsSummaryFile.write(delimiter + delimiter.join(datasetNames) + newline)

        for statisticName in sorted(list(statisticNames)):
            statsSummaryFile.write(statisticName)
            for datasetName in datasetNames:
                statsSummaryFile.write(delimiter + statisticsData[datasetName][statisticName])
            statsSummaryFile.write(newline)

    # Read Library Counts
    guideSequences = set()
    libraryCountData = {}
    for libraryCountFileName in libraryCountFileNames:
        print('Reading library file ' + str(libraryCountFileName))
        key = libraryCountFileName.replace('.library_count.csv','')
        libraryCountData[key] = {}
        with open(join(resultDirectory, libraryCountFileName),'r') as libraryCountFile:
            for line in list(libraryCountFile.readlines()):
                sequence, count = line.strip().split(',')
                libraryCountData[key][sequence] = count.strip()
                guideSequences.add(sequence)

    # Write Library Combined
    libraryCombinedFileName = join(resultDirectory, 'LibraryCountsCombined.csv')
    with open(libraryCombinedFileName,'w') as libraryCombinedFile:
        datasetNames = list(libraryCountData.keys())
        libraryCombinedFile.write(delimiter + delimiter.join(datasetNames) + newline)

        for guideSequence in sorted(list(guideSequences)):
            libraryCombinedFile.write(guideSequence)
            for datasetName in datasetNames:
                libraryCombinedFile.write(delimiter + libraryCountData[datasetName][guideSequence])
            libraryCombinedFile.write(newline)


    # per-gene Statistics
    geneSummaryData = {}
    geneNames = sorted(list(guideRnas.keys()))

    for geneName in geneNames:
        geneSummaryData[geneName] = {}
        for datasetName in datasetNames:
            geneSummaryData[geneName][datasetName] = {}

            geneSummaryData[geneName][datasetName]['read_count'] = 0
            geneSummaryData[geneName][datasetName]['distinct_guides_found'] = 0

            for guideSequence in guideRnas[geneName]:
                geneSummaryData[geneName][datasetName]['read_count'] += int(libraryCountData[datasetName][guideSequence])

                if int(libraryCountData[datasetName][guideSequence]) > 0:
                    geneSummaryData[geneName][datasetName]['distinct_guides_found'] += 1

            if geneSummaryData[geneName][datasetName]['read_count'] == 0:
                print('Warning! for dataset ' + str(datasetName) + ' i did not find any reads for gene ' + str(geneName))

    # Write a per-gene summary
    geneSummaryFileName = join(resultDirectory, 'GeneSummary.csv')
    with open(geneSummaryFileName,'w') as geneSummaryFile:
        for datasetName in datasetNames:
            geneSummaryFile.write(delimiter + datasetName + '_ReadCount')
            geneSummaryFile.write(delimiter + datasetName + '_DistinctGuidesFound')
        geneSummaryFile.write(newline)

        for geneName in geneNames:
            geneSummaryFile.write(geneName)
            for datasetName in datasetNames:
                geneSummaryFile.write(delimiter + str(geneSummaryData[geneName][datasetName]['read_count']))
                geneSummaryFile.write(delimiter + str(geneSummaryData[geneName][datasetName]['distinct_guides_found']))
            geneSummaryFile.write(newline)


if __name__ == '__main__':
    args = parseArgs()

    guideRnas = readGuideRNAs(guideFileName=args.guidernas)

    analyzeDirectory(resultDirectory=args.results, guideRnas=guideRnas)

    print('Done.')