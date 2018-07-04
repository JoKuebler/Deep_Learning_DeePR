# Created by Jonas Kuebler on 05/20/2018
from Bio import SeqIO


# This class takes Fasta File as Input and provides different functionalities such as
# Filtering out sequences according to identifier and write them in to seperate files
# Take any fasta file and write each sequence into seperate file
# Cut out fragments of certain length and write raw output to file
class FastaParser:

    def __init__(self, input_file, output_directory, output_file, identifier=''):
        self.input_file = input_file
        self.identifier = identifier
        self.output_directory = output_directory
        self.output_file = output_file

    # Get fasta sequences with certain identifier by excluding the other sequences
    def filter_by_identifier(self):
        records = list(SeqIO.parse(self.input_file ,'fasta'))
        new_file = open(self.output_directory + self.output_file, 'w')

        for record in records:
            if self.identifier in record.description:
                new_file.write(record)

        new_file.close()

        return new_file

    # Get fasta sequences with certain identifier by excluding the other sequences
    def filter_by_identifier_backwards(self):
        records = list(SeqIO.parse(self.input_file ,'fasta'))
        new_file = open(self.output_directory + self.output_file, 'w')

        for record in records:
            if not any(x in record.description for x in self.identifier):
                new_file.write(record.description + '\n')

        new_file.close()

        return new_file

    # Takes input file and creates file containing fragments of a certain length
    # Excluding certain identifier
    def filter_backwards_cut(self, length, output_file):

        records = list(SeqIO.parse(self.input_file,'fasta'))

        new_file = open(self.output_directory + output_file + '.txt', 'w')

        # Some Peptides hold placeholders for amino acids
        filter_out = ['X', 'x', 'B', 'b', 'Z', 'z', 'U', 'u']

        for record in records:
            if not any(x in record.description for x in self.identifier):
                start = 0
                end = length
                while end < len(record.seq):
                    # Filter out amino acid placeholders
                    if not any(x in record.seq[start:end] for x in filter_out):
                        new_file.write(str(record.seq[start:end]) + '\n')
                    start += length
                    end += length

        new_file.close()

        # Filter out short fragments
        # subprocess.call(["sed", "-r", "-i", '/^.{,33}$/d', str(new_file.name)])

