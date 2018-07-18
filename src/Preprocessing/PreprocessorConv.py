# Created by Jonas Kuebler on 07/15/2018


class PreprocessorConv:

    @staticmethod
    def filter_duplicates(match_file):

        # files with matches
        open_file = open(match_file, 'r')
        
        # list to hold unique matches
        unique_split = []

        # parse each line
        for line in open_file:

            # split at white space
            split = line.split()

            # get rid of path just store filename
            split[1] = split[1].split('/')[-1]

            # replace unnecessary chars
            for ch in ['[', ']', '(', ')']:
                split[2] = split[2].replace(ch, '')
            
            # if filename with identical match positions already exists then filter out
            if (split[1][0:4], split[2]) not in unique_split:
                unique_split.append((split[1][0:4], split[2].replace(']', '')))




