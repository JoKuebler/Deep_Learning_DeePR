# Created by Jonas Kuebler on 05/01/2018
import subprocess

# This class runs TPRpred out of the bioprogs folder on all files given in the input directory
class TPRPredRunner:

    def __init__(self, input_file, output_directory, output_file):
        self.input_file = input_file
        self.output_directory = output_directory
        self.output_file = output_file

    # Runs tpr on fasta file input and stores prediction in output parameters
    def predict_tpr(self):

        run_path = '/ebio/abt1_share/toolkit_sync/bioprogs/tools/tprpred/tprpred'

        # run TPRPred
        subprocess.run([run_path, self.input_file
                        , '-r', '/ebio/abt1_share/toolkit_sync/bioprogs/tools/tprpred/tpr2.8.pp', '-o',
                        str(self.output_directory + self.output_file + '.fa')])