#!/usr/bin/env python

import csv
import sys
import operator as op
import itertools as it

import click

from genome_mapping.validators import Validator
from genome_mapping import mappers
from genome_mapping import matchers


def write(counts):
    def key(overlap):
        return overlap.location.name

    def ordering_key(overlap):
        fn = op.attrgetter('chromosome', 'start', 'stop')
        return fn(overlap.location)

    grouped = it.groupby(sorted(set(counts), key=key), key=key)
    writer = csv.DictWriter(sys.stdout,
                            ['upi',
                             'length',
                             'overlaps',
                             'known_overlaps',
                             'extra_overlaps'])
    writer.writeheader()
    for name, group in grouped:
        overlaps = sorted(group, key=ordering_key)
        if 'URS00002FC3E3' in name:
            from pprint import pprint
            pprint(overlaps, stream=sys.stderr)

        length = (overlaps[0].location.stop - overlaps[0].location.start)
        writer.writerow({
            'upi': name,
            'length': length,
            'overlaps': len(overlaps),
            'known_overlaps': len([g for g in overlaps if g.is_exact()]),
            'extra_overlaps': len([g for g in overlaps if not g.is_exact()]),
        })


@click.command()
@click.argument('genome')
@click.argument('targets')
@click.argument('correct')
@click.option('--matcher', default='exact',
              type=click.Choice(matchers.known()))
@click.option('--mapper', default='blat',
              type=click.Choice(mappers.known()))
def main(genome, targets, correct, matcher='exact', mapper='blat'):
    """
    Try to match some targets to the given genome.
    """

    mapper_class = mappers.fetch(mapper)
    matcher_class = matchers.fetch(matcher)
    validator = Validator(mapper_class(), matcher_class())
    counts = validator(genome, targets, correct)
    write(counts)


if __name__ == '__main__':
    main()
