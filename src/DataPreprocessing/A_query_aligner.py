import subprocess
import os


class Aligner:
    """
    This class is needed for evaluation of the query sequences
    which eventually lead to the positive training data
    """

    def __init__(self, query_dir):

        self.query_dir = query_dir

    def pairwise_align(self):
        """
        Takes directory where query sequences are stored and does
        pairwise alignment
        """

        file_list = os.listdir(self.query_dir)

        # Run TMAlign pairwise with all queries in directory
        for file in file_list:
            for file_two in file_list:
                print(file_two)
                output_file = open('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/queries/TMAlignments/alignment_' + str(file + '_' + file_two) + '.txt', 'w+')
                subprocess.run(['/ebio/abt1_share/update_tprpred/tools/TMalign', '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/queries/third_set/'
                                + file, '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/queries/third_set/' + file_two], stdout=output_file)