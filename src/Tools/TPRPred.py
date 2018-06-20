# Created by Jonas Kuebler on 05/01/2018
import subprocess
import os


# This class runs TPRpred out of the bioprogs folder on all files given in the input directory
class TPRPredRunner:

    def __init__(self, input_file, output_directory):
        self.input_directory = input_file
        self.output_directory = output_directory

    # Runs tpr on each fasta file in input directory
    def predict_tpr(self):

        run_path = '/ebio/abt1_share/toolkit_sync/bioprogs/tools/tprpred/tprpred'

        for file in os.listdir(self.input_directory):

            # run TPRPred
            subprocess.run([run_path, self.input_directory
                            + file.strip(), '-r', '/ebio/abt1_share/toolkit_sync/bioprogs/tools/tprpred/tpr2.8.pp', '-o', str(self.output_directory + str(file[:-3]) + '.fa')])