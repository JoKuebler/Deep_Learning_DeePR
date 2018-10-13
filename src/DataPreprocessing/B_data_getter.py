import os
import re
import warnings
import json
import subprocess
import random
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
        self.trueTprs = {
            '1iyg': [(78, 'RDYVFYLAVGNYRLKEYEKALKYVRGLLQTEPQN')],
            '1hz4': [(14, 'AEFNALRAQVAINDGNPDEAERLAKLALEELPPG'), (53, 'IVATSVLGEVLHCKGELTRSLALMQQTEQMARQH'), (93, 'LWSLIQQSEILFAQGFLQTAWETQEKAFQLINEQ'),
                     (135, 'EFLVRIRAQLLWAWARLDEAEASARSGIEVLSSY'), (174, 'LQCLAMLIQCSLARGDLDNARSQLNRLENLLGNG'), (215, 'SNANKVRVIYWQMTGDKAAAANWLRHTAKPEFAN'),
                     (253, 'QGQWRNIARAQILLGEFEPAEIVLEELNENARSL')],
            '2fez': [(115, 'FVAEKTAGVHAAAAGRFEQASRHLSAALREWRGP'), (171, 'VLAHTAKAEAEIACGRASAVIAELEALTFEHPYR'), (205, 'EPLWTQLITAYYLSDRQSDALGAYRRVKTTLADD')],
            '2aw6': [(191, 'QTVLKNALTISIMNRNLKEAQYYINQFEHLKTIK')],
            '2ond': [(203, 'PEYVLAYIDYLSHLNEDNNTRVLFERVLTSGSLP')],
            '2hr2': [(11, 'AYLALSDAQRQLVAGEYDEAAANCRRAMEISHTM'), (57, 'AFCHAGLAEALAGLRSFDEALHSADKALHYFNRR'), (102, 'ISAVYSRALALDGLGRGAEAMPEFKKVVEMIEER')],
            '4lct': [(108, 'RMGYNDFGDFYYACGMLGDAFKNYIRTRDYCTTT'), (145, 'IHMCMNAILVSIEMGQFTHVTSYVNKAEQNPETL'), (184, 'AKLRCASGLAHLELKKYKLAARKFLDVNPELGNS')]
        }

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
        output_file = open(directory + 'match_dict_updated.json', 'w+')

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
    def write_data(pos_data_arr, filename):

        output_file = open('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/' + filename + '.txt', 'w+')

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
        """
        Filters HHpred results so that only sequences of TPR class are kept
        :param hhr_results: .hhr result files
        :param match_data: match dictionary
        """

        # Dictionary holding all hits for each pdb ID
        new_match_dict = {}

        for file in os.listdir(hhr_results):

            if '.hhr' in file:

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

                        # print(hit)

                        # If hit is found with a prob higher then 50%
                        if self.identifier in hit['hit'] and int(hit['prob'] > 50.0):
                            tpr_fam = True

                # Append correct entry to net dict
                if tpr_fam:

                    subprocess.run(['cp', hhr_results + file, '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/HHpred/results_init_filtered/'])

                    # Get all entries with same chain id
                    for entry in match_data[pdb_id]:
                        if entry["chain"] == chain:

                            # dictionary fill
                            if pdb_id in new_match_dict:
                                new_match_dict[pdb_id].append(entry)
                            else:
                                new_match_dict[pdb_id] = [entry]

                # Remove output file
                # subprocess.run(['rm', hhr_results + 'output.json'])
        # Writes json into query match directory
        self.write_matches_json(self.match_files, new_match_dict)

    def check_range(self, hhr_filtered, match_data):
        """
        Takes alignment and checks if TPR is in same range as true TPR
        :param hhr_filtered: only sequences belonging to TPR class
        :return: confirmed files, unconfirmed files, final raw sequences, final entries as dictionary
        """

        final_seqs = []
        final_entries = {}
        confirmed = []
        unconfirmed = []

        for file in os.listdir(hhr_filtered):

            if '.hhr' in file:

                # Run hhr2json parser and create output file
                output_file = open(hhr_filtered + 'output.json', 'w+')
                subprocess.run(['python3', '/ebio/abt1_share/update_tprpred/code/src/DataPreprocessing/C_hhr2json.py', hhr_filtered + file], stdout=output_file)

                # Open output json
                with open(hhr_filtered + 'output.json') as json_file:

                    # Load data
                    data = json.load(json_file)
                    alignments = data['alignments']

                    # for each file check each alignment
                    for align in alignments:

                        # only consider if alignment hits with higher 50% probability also some alignments are empty
                        if align['info']['probab'] > 50.0 and align['template']['start'] != -1:

                            # query id is the same in all alignments
                            query_id = align['query']['name'].split(':')[0]
                            chain = align['query']['name'][5]

                            # template is the id of the sequence the query was aligned with
                            # one of the seven true tpr sequences
                            template = align['template']['name'][0:4].lower()
                            true_tprs = self.trueTprs[template]
                            # get sequence hit in template
                            templ_seq = align['template']['seq']

                            # for each true tpr
                            for tpr in true_tprs:

                                # check if the tpr is in sequence of template hit
                                if tpr[1] in templ_seq:

                                    # store index
                                    templ_index = templ_seq.index(tpr[1])
                                    # get sequence at this index in query for that hit
                                    query_seq = align['query']['seq'][templ_index:templ_index + 34]

                                    # check if the sequence in the query which aligns to true tpr in template is in positive data
                                    # that means check in match.json
                                    for entry in match_data[query_id]:

                                        # chain has to be the same and if it is the same then take sequence and store it
                                        if entry['chain'] == chain and entry['seq'] == query_seq:

                                            # to not get same fragments check if it is already stored
                                            if query_seq not in final_seqs:

                                                # store confirmed fragment
                                                final_seqs.append(query_seq)

                                                # store confirmed file name
                                                if query_id + '_' + chain not in confirmed:
                                                    confirmed.append(file[:-4])

                                                # dictionary fill to creat updated match.json
                                                if query_id in final_entries:
                                                    final_entries[query_id].append(entry)
                                                else:
                                                    final_entries[query_id] = [entry]

        # Store unconfirmed files
        for file in os.listdir(hhr_filtered):
            if file[:-4] not in confirmed:
                unconfirmed.append(file[:-4])

        return confirmed, unconfirmed, final_seqs, final_entries

    def get_blast_seqs(self, json_dir, match_data):
        """
        After PSI Blast search with all full sequences of initial positive training set this method
        cuts out the sequence in the hit in the range of the true tpr
        :param json_dir: Directory of PSIBlast json result files
        :param match_data: match dictionary
        :return: writes data array to file
        """

        collected_data = []

        # Check all PSI Blast json result files
        for file in os.listdir(json_dir):

            print(file)

            with open(json_dir + file) as json_file:

                # load json and data
                data = json.load(json_file)
                results = data['BlastOutput2'][0]['report']['results']['iterations'][0]['search']
                query_id = results['query_title'][0:4]
                chain = results['query_title'][5]
                hits = results['hits']

                # look at all hits
                for hit in hits:

                    # sequence of the query
                    query_seq = hit['hsps'][0]['qseq']
                    hit_seq = hit['hsps'][0]['hseq']

                    # Check each positive TPR in match_data
                    for entry in match_data[query_id]:

                        # Has to be same chain
                        if entry['chain'] == chain:
                            # tpr of match data
                            tpr = entry['seq']

                            # if tpr is in query seq of hit then take tpr from hit sequence at same position
                            if tpr in query_seq:

                                # get position of tpr in query seq
                                idx = query_seq.index(tpr)

                                # cut out sequence in hit seq on same position
                                hit_tpr = hit_seq[idx:idx + 34]

                                # collect tpr only if the same was not already found and if it has no gaps
                                if hit_tpr not in collected_data and '-' not in hit_tpr:
                                    collected_data.append(hit_tpr)

                print(len(collected_data))
                json_file.close()

            self.write_data(collected_data, 'enriched_pos_data.txt')

    def get_neg_data(self, amount):
        """
        Takes scope70 and cuts all non tpr sequences into 34 long peptides
        :param amount: amount of sequences desired
        :return: writes data array to file
        """

        records = list(SeqIO.parse('/ebio/abt1_share/toolkit_sync/databases/hh-suite/scope70/scope70.fas', 'fasta'))
        collected = []
        x = 0

        rand_sample = random.sample(range(0, len(records)), len(records))

        # go through all records
        for rand in rand_sample:

            if len(collected) > amount:
                break

            # Do not take out of TPR class and sequences smaller then 34 amino acids
            if '118.8' not in records[rand].description and len(records[rand].seq) > 34:

                i = 0

                while i < len(records[rand].seq)-34:

                    sample = records[rand].seq[i:i+34]
                    i += 3

                    if sample not in collected:
                        collected.append(str(sample))

                print(len(collected))

            x += 1

        print(len(collected))
        self.write_data(collected, 'enriched_neg_data')