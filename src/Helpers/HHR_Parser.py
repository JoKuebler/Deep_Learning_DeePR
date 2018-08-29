import os
import subprocess
import json


class HhrParser:

    def __init__(self, hhr_directory=''):
        self.hhr_directory = hhr_directory
        self.identifier = ['a.7.14', 'a.7.16', 'a.23.4', 'a.24.24', 'a.118.5', 'a.118.7', 'a.118.8',
                           'a.118.20', 'a.246.2', 'd.157.1', 'd_26.1' 'e.61.1']

    # keeps only files which are tpr like class in scope
    def filter_files(self, matches_dict):

        counter = 0

        print(len(matches_dict))

        # For all hhr files in directory
        for filename in os.listdir(self.hhr_directory):

            if filename.startswith('hhsearch'):

                name_fasta = filename.split('_')[1]

                if name_fasta in matches_dict:

                    counter = counter + 1
                    print(counter)
                    # store hhr in json file in the directory
                    output_file = open(self.hhr_directory + 'output.json', 'w+')

                    subprocess.run(['python3', '/ebio/abt1_share/update_tprpred/code/src/Helpers/hhr2json.py', self.hhr_directory + filename], stdout=output_file)

                    # open each json file and check if hit is in identifier
                    with open(self.hhr_directory + 'output.json') as file:

                        data = json.load(file)

                        # get PDB Id to find in matches file
                        pdb_id = data['info']['Query'][0:4]
                        chain = data['info']['Query'][5]

                        hits = data['hits']
                        tpr_hit = False

                        # for each hit print the class which is found under 'hit' in json
                        for hit in hits:

                            # set flag to true if hit is in one of the TPR classes
                            if any(x in hit['hit'] for x in self.identifier):

                                template_begin = hit['template_begin']
                                template_end = hit['template_end']

                                if pdb_id in matches_dict:
                                    # for each hit check each entry in match dict and see if hit is in template of hhr
                                    for entry in matches_dict[pdb_id]:

                                        # start pos of match in match dict
                                        start = int(entry['tpr_start'])

                                        # if start pos of match is in template then tpr hit is good to use
                                        if template_begin <= start <= template_end:
                                            tpr_hit = True
                                            continue

                        # if flag is still negative remove the hrr file since the sequence was no TPR hit
                        if not tpr_hit:
                            print('Removed', str(filename))
                            subprocess.run(['rm', self.hhr_directory + filename])

                            # if pdb_id in matches_dict:
                            #
                            #     for index, entry in enumerate(matches_dict[pdb_id]):
                            #         if entry['chain'] == chain:
                            #             del matches_dict[pdb_id][index]
                            #             if len(matches_dict[pdb_id]) == 0:
                            #                 del matches_dict[pdb_id]

                    subprocess.run(['rm', self.hhr_directory + 'output.json'])

                else:
                    subprocess.run(['rm', self.hhr_directory + filename])

        self.write_matches_json(self.hhr_directory, matches_dict)
        print(len(matches_dict))

        return self.hhr_directory

    @staticmethod
    def write_matches_json(directory, match_dict):

        # store hhr in json file in the directory
        output_file = open(directory + 'match_dict.json', 'w+')

        output_file.write(str(match_dict).replace('\'', '"'))

        output_file.close()

    @staticmethod
    def read_matches_json(match_json):

        with open(match_json) as file:
            data = json.load(file)

        return data

    # Looks at hhr directory after filtering and gets fasta sequences -> training data
    def get_fasta_from_hhr(self, fasta_directories):

        fasta_seen = []
        # hhr files
        for filename in os.listdir(self.hhr_directory):

            # get name of pdb id and chain
            pdb_id_chain = filename.split('.')[0][-6:]

            # search for fasta file with same name in directories
            for directory in fasta_directories:

                for fasta_file in os.listdir(directory):

                    if fasta_file[0:6].upper() == pdb_id_chain.upper():

                        if fasta_file not in fasta_seen:

                            fasta_seen.append(fasta_file)

                            # copy to final directory
                            subprocess.run(['cp', directory + fasta_file,
                                            '/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/final_fasta/'])

        return fasta_seen

    @staticmethod
    def proof(matches_dict, final_dir):

        for filename in os.listdir(final_dir):

            pdb_id = filename[0:4]
            print(matches_dict[pdb_id])


if __name__ == '__main__':

    parser = HhrParser('/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/all_hhr/')

    parser.get_fasta_from_hhr(['/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/single_chains_all/'])

    matches_dict = parser.read_matches_json('/ebio/abt1_share/update_tprpred/data'
                                            '/PDB_Approach/All_at_once/tpr_info.json')

    parser.proof(matches_dict, '/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/final_fasta/')