
class PositionParser:

    def __init__(self, alignmentDirectory, goodAlignments):
        self.alignmentDirectory = alignmentDirectory
        self.goodAlignments = goodAlignments

    def findStartingPos(self):

        for file in self.goodAlignments:
            print(file)
            with open(self.alignmentDirectory + file, 'r') as openFile:
                for line in openFile:
                    if line.startswith('("'):
                        print(line)



