# WINDOW CREATER AND ALIGNMENT
# pdbFiles = '/ebio/abt1_share/update_tprpred/data/tpr_like_pdbs/'
# query_sequence = 'data/cAlphas/cAlphas_d1wao11_1TPR.pdb'
# window_directory = '/ebio/abt1_share/update_tprpred/data/windows/'
# alignment_directory = '/ebio/abt1_share/update_tprpred/data/tm_alignments/'
#
# # for each PDB file read it in filter CAlpha and create windows
# for file in os.listdir(pdbFiles):
#
#     # Parse any given PDB File and create File with only C-Alpha atoms
#     PDB_Parser = selfMadeParser(pdbFiles, file)
#
#     # Read in file Content
#     fileContent = PDB_Parser.read_in_structure()
#     print('FILE SUCCESSFULLY READ IN: ', fileContent.name)
#
#     # filter out C-Alpha atoms
#     outputFilePath = PDB_Parser.filterCAlphas(fileContent)
#     print('FILE ONLY CONTAINING CALPHA ATOMS: ', outputFilePath)
#
#     # create window aligner with parameters used for window creation and alignment
#     window_Aligner = windowAligner(outputFilePath, query_sequence, window_directory, alignment_directory)
#
#     # create windows with given length
#     window_Aligner.create_window_files(34)
#
# # align all windows against query sequence
# window_Aligner.align_windows()


# BIOPYTHON FASTA PARSER

# input_fasta = '/ebio/abt1_share/toolkit_sync/databases/hh-suite/scope70/scope70.fas'
# output_directory = '/ebio/abt1_share/update_tprpred/data/tpr_fasta_files/'
# identifier = 'a.118.8'
#
# # create Biopython Fasta Parser
# fasta_Parser = FastaParser(input_fasta, identifier, output_directory)
#
# fasta_Parser.write_to_files()


# Preprocessing methods

# to train neural network we have to encode amino acids
# def encodeAminoAcidSequence(self, fragments):
#     fragments_encoded = []
#
#     for fragment in fragments:
#         joined_string = self.encodeFragment(fragment)
#
#         fragments_encoded.append(joined_string)
#
#
# # takes a protein fragment and encodes it
# @staticmethod
# def encodeFragment(fragment):
#     # Array containing all 20 amino acids
#     amino_acids_array = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C',
#                          'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
#     # empty null vector
#     empty_vector = ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
#                     '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']
#
#     # empty_vector = [False, False, False, False, False, False, False, False, False, False,
#     #                 False,False, False, False, False, False, False, False,False, False]
#
#     # need copy
#     current_vector = empty_vector.copy()
#
#     # array to store vector for each amino acid
#     fragment_vector = []
#
#     # encode each amino acid
#     for aminoAcid in fragment:
#         # find index in amino acid array and set char at index to 1
#         current_vector[amino_acids_array.index(aminoAcid.upper())] = '1'
#
#         # join to get string in the end
#         joined_vector = ','.join(current_vector)
#
#         # append to vector for whole fragment
#         fragment_vector.append(np.array(joined_vector))
#
#         # reset current vector for next amino acied
#         current_vector = empty_vector.copy()
#
#     # join all amino acids to one string
#     joined_string = ':'.join(str(x) for x in fragment_vector)
#
#     return joined_string

# takes a encoded fragment and decodes it
# @staticmethod
# def decodeFragment(encoded_fragment):
#     amino_acids_array = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C',
#                          'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
#
#     # get fragments to decode each aminoacid one by one
#     fragment_list = [encoded_fragment[start:start + 20] for start in range(0, len(encoded_fragment), 20)]
#
#     decoded_fragment = ''
#
#     # decode each encoded amino acid
#     for amino_acid in fragment_list:
#         decoded_fragment += amino_acids_array[amino_acid.index('1')]
#
#     return decoded_fragment