import subprocess
import os


class TMAligner:

    def __init__(self, queryPDB):
        self.queryPDB = queryPDB

    # gets TPR like files from PDB foler in ebio
    # def getTPRFiles(self):
    #     outputFile = open('/ebio/abt1/jkuebler/Desktop/Masterarbeit/TrainingData/TPR_likeIDs.txt', 'r+')
    #     subprocess.run(['grep', '-rl', '/ebio/abt1_share/toolkit_sync/databases/hh-suite/scope70/pdb/', '-e', 'a.118'], stdout=outputFile)
    #
    #     fileList = []
    #
    #     with outputFile as f:
    #         for line in f:
    #             fileList.append(line)
    #
    #     return fileList

    # calls TM aligner for all sequences to align with query sequence
    def alignSequences(self):

        counter = 0

        for file in os.listdir('copiedFiles/'):
            counter = counter + 1
            outputFile = open('/ebio/abt1/jkuebler/Desktop/Masterarbeit/TrainingData/TMAlignments/Alignment_' + str(counter) + '.txt', 'w+')
            subprocess.run(['/ebio/abt1/jkuebler/Desktop/Masterarbeit/TrainingData/TMalign', 'copiedFiles/' + file.strip(), self.queryPDB], stdout=outputFile)

    # gets length of aligned to determine quality
    def getAlignedLength(self, alignmentDirectory):

        finalFileList = []

        for file in os.listdir(alignmentDirectory):
            openFile = open(alignmentDirectory + file, 'rt')

            for line in openFile:
                if line.startswith('TM-score') and 'Chain_2' in line:
                    print(line)
                # content = openFile.read()

                # index = content.find('if normalized by length of Chain_2')
                # print(index)
                # TMScore = content[index:index+20]
                # length = TMScore.split()[len(TMScore.split())-1]
            #
            # if int(length) > 33:
            #     finalFileList.append(file)

        return finalFileList

    # def copyFiles(self, fileList):
    #     print(fileList)
    #
    #     for file in fileList:
    #         subprocess.run(['cp', file.strip(), '/ebio/abt1/jkuebler/Desktop/Masterarbeit/TrainingData/copiedFiles/'])
    #
    #








