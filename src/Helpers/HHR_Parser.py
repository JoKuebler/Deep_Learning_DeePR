import os
import subprocess
import json


class HhrParser:

    def __init__(self, hhr_directory=''):
        self.hhr_directory = hhr_directory
        self.identifier = ['a.7.14', 'a.7.16', 'a.23.4', 'a.24.24', 'a.118.5', 'a.118.7', 'a.118.8',
                           'a.118.20', 'a.246.2', 'd.157.1', 'e.61.1']

    # keeps only files which are tpr like class in scope
    def filter_files(self, matches_dict):

        print(matches_dict)

        # For all hhr files in directory
        for filename in os.listdir(self.hhr_directory):

            # store hhr in json file in the directory
            output_file = open(self.hhr_directory + 'output.json', 'w+')

            subprocess.run(['python3', '/ebio/abt1_share/update_tprpred/code/src/Helpers/hhr2json.py', self.hhr_directory + filename], stdout=output_file)

            # open each json file and check if hit is in identifier
            with open(self.hhr_directory + 'output.json') as file:
                data = json.load(file)

                print(filename)
                # get PDB Id to find in matches file
                pdb_id = data['info']['Query'][0:4]

                hits = data['hits']
                tpr_hit = False

                # for each hit print the class which is found under 'hit' in json
                for hit in hits:

                    # set flag to true if hit is in one of the TPR classes
                    if any(x in hit['hit'] for x in self.identifier):

                        template_begin = hit['template_begin']
                        template_end = hit['template_end']

                        # for each hit check each entry in match dict and see if hit is in template of hhr
                        for entry in matches_dict[pdb_id]:

                            # start pos of match in match dict
                            start = int(entry[0].split(',')[0])

                            # if start pos of match is in template then tpr hit is good to use
                            if template_begin <= start <= template_end:
                                tpr_hit = True
                                continue

                # if flag is still negative remove the hrr file since the sequence was no TPR hit
                if not tpr_hit:
                    print('Removed', str(filename))
                    subprocess.run(['rm', self.hhr_directory + filename])

        subprocess.run(['rm', self.hhr_directory + 'output.json'])


if __name__ == '__main__':

    parser = HhrParser('/ebio/abt1_share/update_tprpred/data/PDB_Approach/Fasta1qqe228/results/hhr_test/')

    parser.filter_files()