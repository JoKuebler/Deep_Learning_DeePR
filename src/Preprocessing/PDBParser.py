# Created by Jonas Kuebler on 05/01/2018
import os

from Bio.PDB import PDBParser, PDBIO


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
    def single_chains(self, directory):

        for file in os.listdir(directory):
            io = PDBIO()
            pdb = PDBParser().get_structure(str(file[0:4]), directory + '/' + file)

            for chain in pdb.get_chains():
                io.set_structure(chain)
                io.save(pdb.get_id() + '_' + chain.get_id() + '.pdb')


