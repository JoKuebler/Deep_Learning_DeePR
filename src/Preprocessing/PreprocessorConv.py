# Created by Jonas Kuebler on 07/15/2018
import subprocess


class PreprocessorConv:

    @staticmethod
    def filter_duplicates(match_file, rmsd_treshold):

        # files with matches
        open_file = open(match_file, 'r')
        
        # list to hold unique matches
        unique_split = []

        unique_dict = {}

        # parse each line
        for line in open_file:

            # split at white space
            split = line.split()

            # get rid of path just store filename
            split[1] = split[1].split('/')[-1]

            # PDB Id in file
            pdb_id = split[1][0:4]

            # RMSD in file
            rmsd = float(split[0])

            # replace unnecessary chars
            for ch in ['[', ']', '(', ')']:
                split[2] = split[2].replace(ch, '')

            tpr_pos = split[2]

            # only keep the ones lower as the treshold
            if rmsd < rmsd_treshold:

                # if filename with identical match positions already exists then filter out
                if (pdb_id, tpr_pos) not in unique_split:

                    unique_split.append((pdb_id, (tpr_pos, rmsd)))

                    # dictionary fill
                    if pdb_id in unique_dict:
                        unique_dict[pdb_id].append((tpr_pos, rmsd))
                    else:
                        unique_dict[pdb_id] = [(tpr_pos, rmsd)]

        return unique_dict

    # Download Fasta Files for PDB IDs
    @staticmethod
    def download_fasta(unique_dict):

        for key in unique_dict:

            subprocess.run(['curl', '-o', '/ebio/abt1_share/update_tprpred/data/PDB_Approach/Fasta/' + str(key) +
                            '.fasta', 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=' +
                            str(key) + '&compressionType=uncompressed'])
