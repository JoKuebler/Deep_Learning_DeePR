import os
import subprocess
import json


class Evaluator:

    def __init__(self):

        self.result_dir = '/ebio/abt1_share/update_tprpred/data/Convolutional/TestFiles/results/'
        self.hhpred_dir = '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/HHpred/yeast_run/e-01/'

    def count_tpr_proteins(self):

        for file in os.listdir(self.result_dir):
            print(file)

            if file.startswith('Pfam'):

                hitlist = []

                with open(self.result_dir + file) as result:

                    for line in result:
                        if line.startswith('HITS'):
                            split = line.split(' ')
                            hits = split[-1].strip()
                            hitlist.append(hits)

                tpr_prots = [hits for hits in hitlist if int(hits) > 1]
                print(len(tpr_prots))

    def count_tprpred_proteins(self):

        for file in os.listdir(self.result_dir):

            if file.startswith('tprpred'):

                print(file)

                tpr_prots = []

                with open(self.result_dir + file) as result:

                    for line in result:

                        if line.startswith('>'):
                            split = line.split('(')
                            if split[0] not in tpr_prots:
                                tpr_prots.append(split[0][1:])

                print(len(tpr_prots))

        return tpr_prots

    def compare(self):

        hit_prot = []

        for file in os.listdir(self.hhpred_dir):

            # Run hhr2json parser and create output file
            output_file = open(self.hhpred_dir + 'output.json', 'w+')
            subprocess.run(['python3', '/ebio/abt1_share/update_tprpred/code/src/DataPreprocessing/C_hhr2json.py', self.hhpred_dir + file], stdout=output_file)

            # Open output json
            with open(self.hhpred_dir + 'output.json') as json_file:

                # Load data
                data = json.load(json_file)

                # get PDB Id to find in matches file
                hits = data['hits']

                for hit in hits:

                    if int(hit['prob'] > 50.0) and hit['struc'] not in hit_prot:

                        hit_prot.append(hit['struc'])

        print(len(hit_prot))
        return hit_prot


if __name__ == '__main__':

    evaluator = Evaluator()

    evaluator.count_tpr_proteins()
    # tpr_found = evaluator.count_tprpred_proteins()
    # hhpred_found = evaluator.compare()
    #
    # print('TPRpred FOUND: ', len(tpr_found))
    # print('HHpred FOUND: ', len(hhpred_found))
    #
    # print(tpr_found)
    #
    # for found in hhpred_found:
    #     if found not in tpr_found:
    #         print(found)