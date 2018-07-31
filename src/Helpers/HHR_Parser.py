import os
import subprocess
import json


class HhrParser:

    def __init__(self, hhr_directory):
        self.hhr_directory = hhr_directory
        self.identifier = ['a.7.14', 'a.7.16', 'a.23.4', 'a.24.24', 'a.118.5', 'a.118.7', 'a.118.8',
                           'a.118.20', 'a.246.2', 'd.157.1', 'e.61.1']

    def filter_files(self):

        # For all hhr files in directory
        for filename in os.listdir(self.hhr_directory):

            # store hhr in json file in the directory
            output_file = open(self.hhr_directory + 'output.json', 'w+')

            subprocess.run(['python3', './hhr2json.py', self.hhr_directory + filename], stdout=output_file)

            # open each json file and check if hit is in identifier
            with open(self.hhr_directory + 'output.json') as file:
                data = json.load(file)

                hits = data['hits']
                tpr_hit = False

                # for each hit print the class which is found under 'hit' in json
                for hit in hits:

                    # set flag to true if hit is in one of the TPR classes
                    if any(x in hit['hit'] for x in self.identifier):
                        tpr_hit = True

                # if flag is still negative remove the hrr file since the sequence was no TPR hit
                if not tpr_hit:
                    subprocess.run(['rm', self.hhr_directory + filename])


if __name__ == '__main__':

    parser = HhrParser('/ebio/abt1_share/update_tprpred/data/PDB_Approach/Fasta1qqe228/results/hhr/')

    parser.filter_files()