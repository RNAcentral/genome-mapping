"""This module contains the Mappers which will map query sequences against
larger target sequence (RNA sequences mapped to genomes or chromosomes). All
Mappers are expected to be callable objects. """

from __future__ import division

import re
import abc
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


def known_formats():
    return {m.format: m for m in ut.children_of(sys.modules[__name__], Mapper)}


def from_format(filename, target_file, format):
    mappers = known_formats()
    mapper = mappers[format]()
    return mapper.parse_result_file(filename, target_file)


class Mapper(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def run(self, genome_file, query_file):
        pass

    @abc.abstractmethod
    def is_valid_sequence(self, sequence):
        return False

    @abc.abstractproperty
    def name(self):
        pass

    @abc.abstractproperty
    def format(self):
        pass

    def valid_sequences(self, query_file, as_dna=False):
        for sequence in SeqIO.parse(query_file, 'fasta'):
            if as_dna:
                sequence.seq = sequence.seq.back_transcribe()

            if self.is_valid_sequence(sequence):
                yield sequence

    def parse_result_file(self, result_file, target_file):
        data = SearchIO.index(result_file, self.format)
        return self.create_mappings(target_file, data)

    def create_mappings(self, query_file, result_index):
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

        for record in self.valid_sequences(query_file):

            urs = re.sub('_\d+$', '', record.id)
            sequence = gm.SequenceSummary(
                id=record.id,
                urs=urs,
                header=record.description,
                length=len(record),
            )

            for hit in result_index[sequence.id]:
                for hsp_index, hsp in enumerate(hit):
                    subhits = []
                    for frag_index, fragment in enumerate(hsp):

                        frag_length = gm.PairStat(
                            query=fragment.query_span,
                            hit=fragment.hit_span,
                        )
                        assert frag_length.hit >= 0, "Bad %s" % frag_length
                        assert frag_length.query >= 0, "Bad %s" % frag_length

                        frag_completeness = gm.PairStat(
                            query=frag_length.query / sequence.length,
                            hit=-1,
                        )

                        start, end = sorted([hsp.hit_start, hsp.hit_end])

                        name = "{urs} ({cur_hsp}/{total_hsp}) ({cur_frag}/{total_frag})".format(
                            urs=sequence.urs,
                            cur_hsp=hsp_index + 1,
                            total_hsp=len(hit),
                            cur_frag=frag_index + 1,
                            total_frag=len(hsp))

                        subhits.append(gm.Fragment(
                            name=name,
                            chromosome=hit.id,
                            start=start,
                            stop=end,
                            is_forward=fragment.hit_strand == 1,
                            stats=gm.FragmentStats(
                                length=frag_length,
                                completeness=frag_completeness,
                            )))

                    assert len({h.is_forward for h in subhits}) == 1
                    assert len(set(h.chromosome for h in subhits)) == 1
                    assert 0 < sum(h.stats.length.query for h in subhits) <= \
                        sequence.length

                    complete = (h.stats.completeness.query for h in subhits)
                    complete = round(sum(complete), 1)
                    assert 0.0 <= complete <= 1.0, "Overly complete: %s" % complete

                    yield gm.Hit(
                        urs=sequence.urs,
                        chromosome=hit.id,
                        start=hsp.hit_start,
                        stop=hsp.hit_end,
                        fragments=subhits,
                        is_forward=subhits[0].is_forward,
                        input_sequence=sequence,
                        stats=gm.Stats.from_summation(
                            [s.stats for s in subhits],
                            hsp.ident_num,
                            gm.PairStat(
                                query=hsp.query_gap_num,
                                hit=hsp.hit_gap_num
                            ),
                        ),
                    )

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
        return self.create_mappings(query_file, output)


class BlatMapper(Mapper):
    """
    This is a simple mapper to map using the BLAT search tool. BLAT is a
    fast method for finding nearly (> 90%) exact matches to a large sequence.
    This uses BioPython to parse the results.
    """

    name = 'blat'
    format = 'blat-psl'

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

    def is_valid_sequence(self, sequence):
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

        options = sorted(set(options + self.default_options))
        with tempfile.NamedTemporaryFile(suffix='.psl') as psl:
            with tempfile.NamedTemporaryFile(suffix='.fa') as query:
                sequences = list(self.sequences(query_path))
                SeqIO.write(sequences, query, 'fasta')
                query.flush()
                cmd = [
                    self.path,
                    '-t=dna',
                    'q=rna',
                ] + options + [
                    genome_file,
                    query.name,
                    psl.name,
                ]
                with open('/dev/null', 'wb') as null:
                    sp.check_call(cmd, stderr=null, stdout=null)
                return SearchIO.index(psl.name, self.format)


class BlastMapper(Mapper):
    name = 'blast'
    format = 'blast-xml'

    default_options = [
        '-word_size=4',
        '-evalue=1',
    ]

    def __init__(self, path='blastn'):
        super(BlastMapper, self).__init__()
        self.path = path

    def is_valid_sequence(self, sequence):
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
        with tempfile.NamedTemporaryFile(suffix='.%s' % self.format) as tmp:
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
                return SearchIO.index(tmp.name, self.format)


class ExonerateMapper(Mapper):
    name = 'exonerate'
    format = 'exonerate-vulgar'

    default_options = [
        # '--model',
        # 'affine:bestfit',
        '--saturatethreshold',
        '0',
    ]

    def __init__(self, path='exonerate'):
        super(ExonerateMapper, self).__init__()
        self.path = path

    def is_valid_sequence(self, sequence):
        return True

    def run(self, genome_file, query_path, options=[]):
        options = self.default_options
        with tempfile.NamedTemporaryFile(suffix='.%s' % self.format) as tmp:
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
                return SearchIO.index(tmp.name, self.format)
