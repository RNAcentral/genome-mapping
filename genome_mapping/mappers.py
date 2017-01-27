"""This module contains the Mappers which will map query sequences against
larger target sequence (RNA sequences mapped to genomes or chromosomes). All
Mappers are expected to be callable objects. """

import os
import sys
import tempfile

from Bio import SeqIO
from Bio import SearchIO

import subprocess as sp

from genome_mapping import data as gm
from genome_mapping import utils as ut

MIN_BLAT_SEQ_LEN = 25
"""The minimum length for sequences to use with BLAT."""


def known():
    return ut.names_of_children(sys.modules[__name__], Mapper)


def fetch(name):
    return ut.get_child(sys.modules[__name__], Mapper, name)


class Mapper(object):
    def create_mappings(self, matches):
        """Create the mappings from the given raw data. The mapping objects in
        genome_mapping.data are clearer (to me at least) so I would rather use
        those sequences instead of the ones provided by BioPython. This
        converts the BioPython results to my results.

        Parameters
        ----------
        matches : list
            A lsit of QueryResult objects to convert.

        Returns
        -------
        results : list
            A list of SequenceResults that represent the matches of the queries
            to the given genome. Note that unlike other parsers if a match is
            in several parts it will be represented as several matches in this
            list, that is to say there will be several MappingHit's for it.
        """

        results = {}
        for result in matches:
            if result.id not in results:
                results[result.id] = gm.SequenceResults(name=result.id)
            current = results[result.id]
            for hit in result:
                for fragment in hit:
                    stats = gm.HitStats(identical=fragment.ident_num,
                                        gaps=fragment.gapopen_num,
                                        query_length=result.seq_len,
                                        hit_length=fragment.hit_end - fragment.hit_start,
                                        )
                    current.add(
                        gm.MappingHit(
                            name=result.id,
                            chromosome=hit.id,
                            start=fragment.hit_start,
                            stop=fragment.hit_end,
                            is_forward=fragment[0].hit_strand == 1,
                            stats=stats))
        return results.values()

    def __call__(self, genome_file, query_file):
        """Perform the mapping. This takes a genome_file and a list of
        sequences and produces the mappings for those sequences. Note that if
        there are no valid sequences in the list then None is returned.

        Parameters
        ----------
        genome_file : str
            Path to the genome file.

        sequences : list
            List of all sequences to map. Each object in the list should be
            writable by BioPython.

        Returns
        -------
        mappings : list
            A list of mapping objects.
        """

        # filtered = self.filter_sequences(sequences)
        # if not filtered:
        #     return None
        # query_path = self.create_query(filtered)
        output = self.run(genome_file, query_file)
        return self.create_mappings(output)


class BlatMapper(Mapper):
    """
    This is a simple mapper to map using the BLAT search tool. BLAT is a
    fast method for finding nearly (> 90%) exact matches to a large sequence.
    This uses BioPython to parse the results.
    """

    name = 'blat'

    def __init__(self, path='blat'):
        """Create a new BlatMapper.

        Parameters
        ----------
        path : str
            The full path to the BLAT binary.
        """
        self.path = path

    def filter_sequences(self, sequences):
        """
        Filter all sequences to only those that are long enough. BLAT does
        not work with short,  25 nt, sequences.

        Parameters
        ----------
        sequences : list
            A list of sequences to filter

        Returns
        -------
        valid_sequences : list
            The list of long enough sequences.
        """
        for sequence in sequences:
            if len(sequence) > MIN_BLAT_SEQ_LEN:
                yield sequence

    def create_query(self, sequences):
        """This creates a query file for BLAT. It will write out all given
        sequences to a FASTA file and return the path to that file. This file
        will be used by BLAT to query.

        Parameters
        ----------
        sequences : list
            List of sequences to write.

        Returns
        -------
        path : str
            The path to where the sequences were written.
        """

        filename = os.path.abspath('query-sequences.fa')
        valid = self.filter_sequences(sequences)
        with open(filename, 'wb') as handle:
            SeqIO.write(valid, handle, 'fasta')
        return filename

    def run(self, genome_file, query_path):
        """Run the BLAT program on the given genome with the given query.

        Parameters
        ----------
        genome_file : str
            Full path to the genome file to use.

        query_path : str
            Full path to the query file to use.

        Returns
        -------
        results : list
            A list of parsed QueryResult objects. The parsing is done by
            BioPython.
        """

        with tempfile.NamedTemporaryFile(suffix='.psl') as tmp:
            sp.check_call([self.path, genome_file, query_path, tmp.name])
            return list(SearchIO.parse(tmp.name, 'blat-psl'))
