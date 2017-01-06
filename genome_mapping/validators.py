"""This module contains useful classes for detecting if matches of some
RNA sequences are correct. See Validator for details.
"""

import copy

from genome_mapping import data as gm


class Validator(object):
    """A Validator is a class that is used to check how well a particular
    alignment tool (Mapper) and a particular set of cutoffs for 'good' matches
    (MappingFilter) do on a given data set. It is used to validate if our logic
    for mapping RNA sequences to a genome is good. Standard usage of this is to
    use a genome with some known mappings and attempt to reproduce them. If
    they can be reproduced then the setup used is at least reasonable for that
    genome and set of RNA's.

    This object is highly configurable because the requirements for each genome
    and type of RNA are likely to be different.
    """

    def __init__(self, mapper, matcher):
        """Create a new Validator.

        Parameters
        ----------
        mapper : Mapper
            An object to call that will map FASTA sequences to a target
            sequence.
        matcher : MappingFilter
            A Matcher that will filter out matches from mapping that are not
            valid. For example using an ExactMappingFilter will require that
            all matches between query and target sequences are exact.
        """

        self.mapper = mapper
        self.matcher = matcher

    def __call__(self, genome_file, target_file, given_expected):
        """This will use the objects mapper to map all sequences in the
        target_file to the genome in the genome_file. It will then compare the
        mappings to the mappings in given_expected and produce a summary. Only
        the matches which are valid according to this objects matcher will be
        used.

        Parameters
        ----------
        genome_file : str
            Full path to the file containing the genome to query agains. Must
            be a file that is readable by BLAT.

        target_file : str
            Full path to the file containing the RNA sequences to query with.
            Must be a FASTA file.

        given_expected : dict
            This must be a dict which maps from accession number to a set of
            Mapping objects. Each Mapping represents the correct location in
            the genome for the given sequence.

        Returns
        -------
        counts : SummaryCounts
            A set of counts summarizing the overall matches.
        """

        counts = gm.SummaryCounts()
        mappings = self.mapper(genome_file, target_file)
        if mappings is None:
            raise ValueError("No mappings produced")
        valid = self.matcher.filter_matches(mappings)

        expected = copy.deepcopy(given_expected)
        for mapping in valid:
            key = mapping.accession
            if key not in expected:
                counts.extra_matches += 1
            result = mapping.as_mapping()
            if key in expected:
                correct = expected.pop(key)
                if correct == result:
                    counts.valid_matches += len(correct)
        return counts
