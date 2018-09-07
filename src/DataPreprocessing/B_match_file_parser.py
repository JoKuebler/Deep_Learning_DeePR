import os
import re
from Bio import SeqIO
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


class MatchParser:
    """
    This class parses the match files obtained from running
    queries against PDB with MASTER
    """

    def __init__(self, match_files, all_matches):

        self.match_files = match_files
        self.all_matches = all_matches

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

        # Do for all query match files
        for file in match_files:

            content = open(self.match_files + file)

            # Extracts first line (best match with itself) and gets PDB ID of query
            pdb_id = content.readline().split()[1].split('/')[-1][0:4].upper()

            # Get match files for current query in list
            current_matches = [x for x in os.listdir(self.all_matches) if pdb_id.lower() in x]

            # Iterate through query match file to make sure for every line there is a separate match.pdb file
            for idx, val in enumerate(content):

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
                    match_entry = [{"rmsd": rmsd, "chain": chain_id, "tpr_start": start_pos, "seq": str(seq)}]
                    # dictionary fill
                    if pdb_id_match in match_dict:
                        match_dict[pdb_id_match].append(match_entry[0])
                    else:
                        match_dict[pdb_id_match] = match_entry
                else:
                    wrong.append(str(seq))

        # Writes json into query match directory
        self.write_matches_json(self.match_files, match_dict)

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
        output_file = open(directory + 'match_dict.json', 'w+')

        output_file.write(str(match_dict).replace('\'', '"'))

        output_file.close()
