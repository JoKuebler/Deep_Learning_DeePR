# Created by Jonas Kuebler on 05/20/2018
from Bio import SeqIO


# This class takes Fasta File as Input and provides different functionalities such as
# Filtering out sequences according to identifier and write them in to seperate files
# Take any fasta file and write each sequence into seperate file
# Cut out fragments of certain length and write raw output to file
class FastaParser:

    def __init__(self, input_file, identifier, output_directory):
        self.input_file = input_file
        self.identifier = identifier
        self.output_directory = output_directory

    # writes each record to seperate file
    def write_to_files_by_identifier(self):
        records = list(SeqIO.parse(self.input_file ,'fasta'))

        for record in records:

            if self.identifier in record.description:
                SeqIO.write(record, self.output_directory + record.id + '.fa', 'fasta')

    def write_to_files_all(self):

        records = list(SeqIO.parse(self.input_file, 'fasta'))

        for record in records:
            SeqIO.write(record, self.output_directory + record.id + '.fa', 'fasta')

    # cut out 34 fragments of non tpr like
    def write_to_files_fragments(self, length):

        records = list(SeqIO.parse(self.input_file,'fasta'))
        new_file = open(self.output_directory + 'raw_nontpr' + '.txt', 'w')

        for record in records:
            if self.identifier not in record.description:
                new_file.write(str(record.seq[0:length]) + '\n')
                new_file.write(str(record.seq[34:length+34]) + '\n')
                new_file.write(str(record.seq[68:length+68]) + '\n')
                new_file.write(str(record.seq[102:length+102]) + '\n')
                new_file.write(str(record.seq[136:length + 136]) + '\n')
                new_file.write(str(record.seq[170:length + 170]) + '\n')

        new_file.close()

