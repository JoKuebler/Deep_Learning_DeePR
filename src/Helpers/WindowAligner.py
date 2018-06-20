# Created by Jonas Kuebler on 05/04/2018
import os
import subprocess

global_counter = 0


# This class provides a sliding window alogrithm in order to cut pdb files in windows of desired length
class WindowAligner:

    def __init__(self, calpha_struc_path, query_struc, window_directory, alignment_directory):
        self.cAlphaStructurePath = calpha_struc_path
        self.structure = calpha_struc_path.split('/')[len(calpha_struc_path.split('/')) - 1]
        self.query_struc = query_struc
        self.window_directory = window_directory
        self.alignment_directory = alignment_directory

    # creates file for each window in sliding window approach
    def create_window_files(self, window_size):

        with open(self.cAlphaStructurePath) as fileContent:

            # declare start of window
            window_start = int()

            # declare index of filename
            file_name_index = 1

            # find first line of structure
            for num, line in enumerate(fileContent):
                # skip all the headers and remember position of first structure residue
                if not line.startswith('REMARK'):
                    window_start = num
                    break

            # reset pointer of file
            fileContent.seek(0)

            # read in all lines of file
            read_lines = fileContent.readlines()

            # get length of file in order to access via index later
            file_length = len(read_lines)

            # calculate end of first window by adding size to start position
            window_end = window_start + window_size

            # slide window over file until end is reached
            while window_end <= file_length:

                # create file for each iteration with specific filename index
                created_window_file = open(self.window_directory + self.structure[:-4] + '_window' +
                                           str(file_name_index) + '.pdb', 'w')

                # write lines of sliding window to file
                created_window_file.writelines(read_lines[window_start:window_end])

                # increment file name
                file_name_index += 1

                # increment window start pos
                window_start += 1

                # increment window end pos
                window_end += 1

    # aligns windows to query sequence using TM-align
    def align_windows(self):

        file_index = 1

        # align each file in the window directory
        for file in os.listdir(self.window_directory):

            # store alignment in txt file in the alignment directory
            output_file = open(self.alignment_directory + self.structure[:-4] +
                               '_alignment_window_' + str(file_index) + '.txt', 'w+')

            # run TM aligner
            subprocess.run(['/ebio/abt1/jkuebler/Desktop/Masterarbeit/TrainingData/TMalign', self.window_directory
                            + file.strip(), self.query_struc], stdout=output_file)

            output_file.close()

            # look at each file ...
            open_output_file = open(self.alignment_directory + self.structure[:-4] +
                                    '_alignment_window_' + str(file_index) + '.txt', 'rt')

            # get file content
            file_content = open_output_file.read()

            # find RMSD value
            index = file_content.find('RMSD')
            rmsd = file_content[index+8:index+12]

            # ...and only keep it if RMSD threshold is not exceeded
            if float(rmsd) >= 2:
                os.remove(self.alignment_directory + self.structure[:-4] +
                          '_alignment_window_' + str(file_index) + '.txt')
                os.remove(self.window_directory + file.strip())

            file_index += 1








