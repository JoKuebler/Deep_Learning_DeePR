import os
import re
import warnings
import json
import subprocess
from Bio import SeqIO
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


class DataGetter:
    """
    This class parses the match files obtained from running
    queries against PDB with MASTER
    """

    def __init__(self, match_files, all_matches):

        self.match_files = match_files
        self.all_matches = all_matches
        self.identifier = 'a.118.8'

    def parse_info(self):
        """
        Parse relevant information to obtain sequence for training
        data
        """
        # Match files listing all hits with query
        match_files = os.listdir(self.match_files)
        # Raw sequences of hits
        match_seqs = []
        # Dictionary holding all hits for each pdb ID
        match_dict = {}

        wrong = []

        counter = 1

        # Do for all query match files
        for file in match_files:

            print(file)

            # Consider only match files in folder
            if ".match" in file:

                content = open(self.match_files + file)

                # Extracts first line (best match with itself) and gets PDB ID of query
                pdb_id = content.readline().split()[1].split('/')[-1][0:4].upper()

                # Get match files for current query in list
                current_matches = [x for x in os.listdir(self.all_matches) if pdb_id.lower() in x]

                # Iterate through query match file to make sure for every line there is a separate match.pdb file
                for idx, val in enumerate(content):

                    print(counter)
                    counter += 1
                    # content of current match.pdb file to get position
                    match_content = open(self.all_matches + current_matches[idx]).readlines()
                    # split first and second line
                    first_line_split, second_line_split = match_content[0].split(), match_content[1].split()

                    # RMSD always in first line in second column
                    rmsd = first_line_split[1]
                    # start position is in second line always
                    start_pos = second_line_split[5]
                    # PDB ID of match
                    pdb_id_match = first_line_split[2].split('/')[-1][0:4].upper()

                    # if no integer is found then start position is in different column (e.g. A1167)
                    # then also the chain ID is found in that column instead of the previous one
                    if re.match(r"^\d+$", start_pos):
                        chain_id = second_line_split[4]
                    else:
                        start_pos, chain_id = second_line_split[4][1:], second_line_split[4][0]

                    # Get sequence from match.pdb file and append to all sequences
                    seq = self.get_seq(self.all_matches + current_matches[idx])

                    # If hit is somehow not 34 aa long dont use it
                    # Also makes sure that exact same sequence is not already seen
                    if len(str(seq)) == 34 and str(seq) not in match_seqs:
                        match_seqs.append(str(seq))

                        # Build match dictionary
                        match_entry = [{"rmsd": rmsd, "chain": chain_id, "tpr_start": start_pos, "seq": str(seq), "query": pdb_id}]
                        # dictionary fill
                        if pdb_id_match in match_dict:
                            match_dict[pdb_id_match].append(match_entry[0])
                        else:
                            match_dict[pdb_id_match] = match_entry
                    else:
                        wrong.append(str(seq))

        # Writes json into query match directory
        self.write_matches_json(self.match_files, match_dict)

        print(wrong)

    @staticmethod
    def get_seq(pdb_file):
        """
        Gets sequence from PDB file
        :param pdb_file: pdb file of the match
        :return:
        """
        handle = open(pdb_file, "rU")
        records = list(SeqIO.parse(handle, "pdb-atom"))

        return records[0].seq

    @staticmethod
    def write_matches_json(directory, match_dict):
        """
        Writes dictionary to json
        :param directory:
        :param match_dict: dictionary containing hits
        :return:
        """

        # store hhr in json file in the directory
        output_file = open(directory + 'match_dict_queryID.json', 'w+')

        output_file.write(str(match_dict).replace('\'', '"'))

        output_file.close()

    @staticmethod
    def write_to_fasta(match_data):
        """
        Write each sequence in the match json into fasta format in one fasta file
        Important e.g TPRpred
        :param match_data: dictionary containing hits
        :return:
        """

        output_file = open('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/fasta.fa', 'w+')

        for key in match_data:
            for entry in match_data[key]:

                header = str(">" + key + entry["chain"] + '\n')

                output_file.write(header)
                output_file.write(entry["seq"] + '\n')

        output_file.close()

    @staticmethod
    def read_match_json(match_json):

        pos_data = []

        with open(match_json) as file:
            data = json.load(file)

        for key in data:
            for entry in data[key]:
                pos_data.append(entry["seq"])

        return data, pos_data

    @staticmethod
    def write_pos_data(pos_data_arr):

        output_file = open('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/new_pos_data.txt', 'w+')

        for seq in pos_data_arr:
            output_file.write(seq + '\n')

    @staticmethod
    def download_fasta(match_data, download_dir):
        """
        Takes match dictionary and downloads all fasta files containing hits
        :param match_data:
        :param download_dir:
        :return:
        """

        print('FASTA Files to be downloaded: ' + str(len(match_data)))
        for key in match_data:

            subprocess.run(['curl', '-o', download_dir + str(key) +
                            '.fasta', 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=' +
                            str(key) + '&compressionType=uncompressed'])

    @staticmethod
    def single_chains_fasta(match_data, download_directory, output_directory):
        """
        Takes full length sequences of PDB ID and only keeps the chain a hit is found in
        :param match_data
        :param download_directory:
        :param output_directory:
        :return:
        """

        final_records = []

        for file in os.listdir(download_directory):

            records = list(SeqIO.parse(download_directory + file, 'fasta'))

            for rec in records:
                pdb_id = rec.id[0:4]
                chain = rec.id[5]

                for entry in match_data[pdb_id]:
                    if chain == entry["chain"]:
                        final_records.append(rec)
                        break

        for record in final_records:

            SeqIO.write(record, output_directory + record.name[0:6].replace(':', '_') + '.fasta', "fasta")

    def hhpred_init_filter(self, hhr_results, match_data):

        # Dictionary holding all hits for each pdb ID
        new_match_dict = {}

        for file in os.listdir(hhr_results):

            # Run hhr2json parser and create output file
            output_file = open(hhr_results + 'output.json', 'w+')
            subprocess.run(['python3', '/ebio/abt1_share/update_tprpred/code/src/DataPreprocessing/C_hhr2json.py', hhr_results + file], stdout=output_file)

            # Open output json
            with open(hhr_results + 'output.json') as json_file:

                # for each output file tpr_fam gets reset
                tpr_fam = False

                # Load data
                data = json.load(json_file)

                # get PDB Id to find in matches file
                pdb_id = data['info']['Query'][0:4]
                chain = data['info']['Query'][5]
                # get all hits
                hits = data['hits']

                # Check all hits
                for hit in hits:

                    # If hit is found with a prob higher then 50%
                    if self.identifier in hit['hit'] and int(hit['prob'] > 50.0):
                            tpr_fam = True

            # Append correct entry to net dict
            if tpr_fam:

                # Get all entries with same chain id
                for entry in match_data[pdb_id]:
                    if entry["chain"] == chain:

                        # dictionary fill
                        if pdb_id in new_match_dict:
                            new_match_dict[pdb_id].append(entry)
                        else:
                            new_match_dict[pdb_id] = [entry]

            # Remove output file
            subprocess.run(['rm', hhr_results + 'output.json'])

            print(new_match_dict)