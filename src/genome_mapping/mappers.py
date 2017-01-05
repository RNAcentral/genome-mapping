"""This module contains the Mappers which will map query sequences against
larger target sequence (RNA sequences mapped to genomes or chromosomes). All
Mappers are expected to be callable objects. """

import os

from Bio import SeqIO
from Bio import SearchIO

import subprocess as sp

from genome_mapping import data as gm

MIN_BLAT_SEQ_LEN = 25
"""The minimum length for sequences to use with BLAT."""


class BlatMapper(object):
    """This is a simple mapper to map using the BLAT search tool. BLAT is a
    fast method for finding nearly (> 90%) exact matches to a large sequence.
    This uses BioPython to parse the results.
    """

    def __init__(self, path):
        """Create a new BlatMapper.

        Parameters
        ----------
        path : str
            The full path to the BLAT binary.
        """
        self.path = path

    def filter_sequences(self, sequences):
        """Filter all sequences to only those that are long enough. BLAT does
        not work with short (< 25nt) sequences.

        Parameters
        ----------
        sequences : list
            A list of sequences to filter

        Returns
        -------
        valid_sequences : list
            The list of long enough sequences.
        """
        return [seq for seq in sequences if len(seq) > MIN_BLAT_SEQ_LEN]

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
        with open(filename, 'wb') as handle:
            SeqIO.write(sequences, handle, 'fasta')
        return filename

    def run_blat(self, genome_file, query_path):
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
        output = os.path.abspath('output.psl')
        sp.check_call([self.path, '-q', 'rna', genome_file, query_path,
                       output])
        return SearchIO.read(output, 'blat-psl')

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
            for hit in result:
                for fragment in hit:
                    stats = gm.HitStats(identical=hit.ident_num,
                                        gaps=hit.gapopen_num,
                                        query_length=results.seq_len,
                                        hit_length=hit.seq_len,
                                        )
                    results[result.id].add(
                        gm.MappingHit(
                            chromosome=hit.id,
                            start=hit.hit_start,
                            stop=hit.hit_stop,
                            is_forward=fragment.hit_strand,
                            stats=stats))
        return results

    def __call__(self, genome_file, sequences):
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

        filtered = self.filter_sequences(sequences)
        if not filtered:
            return None
        query_path = self.create_query(filtered)
        output = self.run_blat(query_path)
        return self.create_mappings(output)
