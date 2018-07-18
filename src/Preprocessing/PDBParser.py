# Created by Jonas Kuebler on 05/01/2018

import subprocess
import warnings
from pathlib import Path

from Bio import BiopythonWarning
from Bio.PDB import PDBList
from Bio.PDB import PDBParser, PDBIO

warnings.simplefilter('ignore', BiopythonWarning)

# This class takes an pdb File as input and returns
# a PDB file containing only the CAlpha atoms of the input file
class PDB_Parser:

    # initialize obbject with path to structure, structure name and new file to write in calphas only
    def __init__(self, files_path='', structure_name=''):
        self.files_path = files_path
        self.structure_name = structure_name
        # self.createdFile = open('/ebio/abt1_share/update_tprpred/data/calphas/calphas_' + self.structure_name, 'w')

    # reads in structure
    def read_in_structure(self):

        # open the file
        file_content = open(self.files_path + self.structure_name, 'rt')

        # return openFile
        return file_content

    # filters out C-Alpha atoms because TM align only looks at those
    # def filterCAlphas(self, file_content):
    #
    #     for line in file_content:
    #         if line.startswith('REMARK'):
    #             self.createdFile.write(line)
    #         elif 'CA' in line:
    #             self.createdFile.write(line)
    #
    #     self.createdFile.close()
    #
    #     return self.createdFile.name

    # Takes folder of pdb files and writes each chain into seperate file
    @staticmethod
    def single_chains(file, directory, entry):

        io = PDBIO()
        pdb = PDBParser().get_structure(str(entry), file)

        chain_files = []

        # Check how many chains there are
        for element in pdb.get_chains():
            chain_files.append(str(pdb.get_id()) + '_' + str(element.get_id()))

        # Filter duplicates
        chain_files = list(set(chain_files))

        # Write them to seperate files
        if len(chain_files) > 1:
            for chain in pdb.get_chains():
                io.set_structure(chain)
                io.save(directory + '/' + pdb.get_id() + '_' + chain.get_id() + '.pdb')

        return chain_files

    # Run the createPDS from master
    @staticmethod
    def build_master_database(directory, pdb_file):

        subprocess.run(['/tmp/jonas/createPDS', '--type', 'target', '--pdb',
                        directory + pdb_file + '.pdb', '--pds',
                        directory + pdb_file + '.pds'])

        subprocess.run(['rm', directory + pdb_file + '.pdb'])

    # Downloads all PDB Files and creates Database for master search
    def download_build(self, directory):

        # Init PDBList from Biopython
        pdbl = PDBList()
        all_pdb = pdbl.get_all_entries()
        not_found = []

        with open('/tmp/jonas/pdb_exclude.txt') as file_content:
            leave_out = file_content.readlines()

        leave_out_strip = [x.strip() for x in leave_out]

        # testIDs = ['3hhb', '2hr2', '2hhb', '1zb1', '5MQX','111f', '1y8m', '2pqn', '2pqr', '1nzn', '1pc2']

        # Get all entries
        for entry in all_pdb[11194:]:

            if entry not in leave_out_strip:

                # Download via Biopython API
                pdbl.retrieve_pdb_file(entry, obsolete=False,
                                       pdir=directory, file_format='pdb')

                file = Path(directory + 'pdb' + str(entry).lower() + '.ent')

                if file.is_file():

                    # Rename .ent files to .pdb files
                    subprocess.run(['mv', directory + 'pdb' + str(entry.lower()) + '.ent',
                                    directory + str(entry).lower() + '.pdb'])

                    # Check if pdb file has multiple chains, if yes write to multiple files
                    chains = list(set(self.single_chains(directory + str(entry.lower()) + '.pdb', directory, entry)))

                    # If there are multiple chain files create pds file for all
                    if len(chains) > 1:

                        subprocess.run(['rm', directory + str(entry).lower() + '.pdb'])

                        for chain in chains:
                            self.build_master_database(directory, chain)

                    # pds file for single
                    else:
                        self.build_master_database(directory, str(entry).lower())
                else:
                    print('NOT FOUND')
                    not_found.append(file)
                    continue

    # runs master with parameters specified
    @staticmethod
    def run_master(directory, query, database_list, rmsd):

        subprocess.run(['/tmp/jonas/run/master', '--query', str(directory + query),
                        '--targetList', str(directory + database_list),
                        '--rmsdCut', str(rmsd), '--matchOut', str(directory + 'query.match'),
                        '--seqOut', str(directory + 'query.seq'),
                        '--structOut', str(directory + 'query.struc')])


# Main Method
if __name__ == '__main__':

    # parser_object = PDB_Parser()
    #
    # # parser_object.download_build('/ebio/abt1_share/update_tprpred/data/PDB_Approach/PDS/')
    #
    # parser_object.run_master('/tmp/jonas/run/',
    #                          'query.pds',
    #                          'database_list', 4.0)

    pdbl = PDBList()
    all_pdb = pdbl.get_all_entries()

    # 12:30 9.6GB 15891
    # 13:00 11.GB 16743
    # 15:30 13GB  20908
    # 17:00 15GB  23668
    print(all_pdb.index('3Q33'))
