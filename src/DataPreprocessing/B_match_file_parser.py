import subprocess
import os


class MatchParser:
    """
    This class parses the match files obtained from running
    queries against PDB with MASTER
    """

    def __init__(self, match_dir):

        self.match_dir = match_dir

    def parse_info(self):
        """
        Parse relevant information to obtain sequence for training
        data
        """

        # TODO WRITE PARSER AND USE PDB MATCH FILES