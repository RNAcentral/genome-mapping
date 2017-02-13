"""This module contains the Mappers which will map query sequences against
larger target sequence (RNA sequences mapped to genomes or chromosomes). All
Mappers are expected to be callable objects. """

import re
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
    def __init__(self):
        self.mapping = {}

    def valid_sequence(self, sequence):
        return False

    def lookup_sequence(self, result):
        upi = re.sub('_\d+$', '', result.id)
        return self.mapping[upi]

    def sequences(self, query_file, as_dna=False):
        for sequence in SeqIO.parse(query_file, 'fasta'):
            if as_dna:
                sequence.seq = sequence.seq.back_transcribe()

            if self.valid_sequence(sequence):
                uri = re.sub('_\d+$', '', sequence.id)
                self.mapping[uri] = gm.SequenceSummary(
                    id=sequence.id,
                    uri=uri,
                    header=sequence.description,
                )
                yield sequence

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

        for result in matches:
            sequence = self.lookup_sequence(result)
            for hit in result:
                for fragment in hit:
                    hit_len = fragment.hit_end - fragment.hit_start
                    identity = getattr(fragment, 'ident_pct', None)
                    gaps = getattr(fragment, 'gapopen_num', 0)
                    stats = gm.Stats(identical=fragment.ident_num,
                                     identity=identity,
                                     gaps=gaps,
                                     query_length=result.seq_len,
                                     hit_length=hit_len)

                    start, end = sorted([fragment.hit_start, fragment.hit_end])
                    yield gm.Hit(
                        name=result.id,
                        chromosome=hit.id,
                        start=start,
                        stop=end,
                        is_forward=fragment[0].hit_strand == 1,
                        input_sequence=sequence,
                        stats=stats)

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

        output = self.run(genome_file, query_file)
        return self.create_mappings(output)


class BlatMapper(Mapper):
    """
    This is a simple mapper to map using the BLAT search tool. BLAT is a
    fast method for finding nearly (> 90%) exact matches to a large sequence.
    This uses BioPython to parse the results.
    """

    name = 'blat'

    default_options = [
        '-noTrimA',
        '-fine',
        '-minIdentity=0',
        '-minScore=0',
        '-minMatch=1',
        '-maxGap=3',
    ]

    def __init__(self, path='blat'):
        """Create a new BlatMapper.

        Parameters
        ----------
        path : str
            The full path to the BLAT binary.
        """
        super(BlatMapper, self).__init__()
        self.path = path

    def valid_sequence(self, sequence):
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
        return len(sequence) > MIN_BLAT_SEQ_LEN

    def run(self, genome_file, query_path, options=[]):
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

        options = self.default_options
        with tempfile.NamedTemporaryFile(suffix='.psl') as tmp:
            with tempfile.NamedTemporaryFile(suffix='.fa') as qtmp:
                SeqIO.write(self.sequences(query_path), qtmp, 'fasta')
                cmd = [
                    self.path,
                    '-t=dna',
                    'q=rna',
                ] + options + [
                    genome_file,
                    qtmp.name,
                    tmp.name,
                ]
                with open('/dev/null', 'wb') as null:
                    sp.check_call(cmd, stdout=null)
                return list(SearchIO.parse(tmp.name, 'blat-psl'))


class BlastMapper(Mapper):
    name = 'blast'

    default_options = [
        '-word_size=4',
        '-evalue=1',
    ]

    def __init__(self, path='blastn'):
        super(BlastMapper, self).__init__()
        self.path = path

    def valid_sequence(self, sequence):
        return True

    def run(self, genome_file, query_path, options=[]):
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

        options = self.default_options
        format = 'blast-xml'
        with tempfile.NamedTemporaryFile(suffix='.%s' % format) as tmp:
            with tempfile.NamedTemporaryFile(suffix='.fasta', mode='wb') as qtmp:
                SeqIO.write(self.sequences(query_path), qtmp, 'fasta')
                cmd = [
                    self.path,
                ] + options + [
                    '-outfmt=5',
                    '-db=%s' % genome_file,
                    '-query=%s' % qtmp.name,
                    '-out=%s' % tmp.name,
                ]
                with open('/dev/null', 'wb') as null:
                    sp.check_call(cmd, stdout=null)
                return list(SearchIO.parse(tmp.name, format))


class ExonerateMapper(Mapper):
    name = 'exonerate'

    default_options = [
        # '--model',
        # 'affine:bestfit',
        '--saturatethreshold',
        '0',
    ]

    def __init__(self, path='exonerate'):
        super(ExonerateMapper, self).__init__()
        self.path = path

    def valid_sequence(self, sequence):
        return True

    def run(self, genome_file, query_path, options=[]):
        options = self.default_options
        format = 'exonerate-vulgar'
        with tempfile.NamedTemporaryFile(suffix='.%s' % format) as tmp:
            with tempfile.NamedTemporaryFile(suffix='.fasta') as qtmp:
                SeqIO.write(self.sequences(query_path, as_dna=True), qtmp, 'fasta')
                cmd = [
                    self.path,
                ] + options + [
                    '--showvulgar',
                    'TRUE',
                    # '-E',
                    # 'TRUE',
                    '--querytype',
                    'dna',
                    '--targettype',
                    'dna',
                    '--query',
                    qtmp.name,
                    '--target',
                    genome_file,
                ]
                sp.check_call(cmd, stdout=tmp)
                return list(SearchIO.parse(tmp.name, format))
