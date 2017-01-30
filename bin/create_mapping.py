#!/usr/bin/env python

import itertools as it

import click

from genome_mapping.validators import Validator
from genome_mapping import mappers
from genome_mapping import matchers


def write(counts):
    from pprint import pprint
    key = lambda o: o.location.name
    grouped = it.groupby(sorted(counts, key=key), key=key)
    for name, group in grouped:
        overlaps = list(group)
        exact = any(g.is_exact() for g in overlaps)
        if not exact:
            pprint(overlaps)


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
