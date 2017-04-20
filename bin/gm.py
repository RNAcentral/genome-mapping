#!/usr/bin/env python -W ignore

import os
import csv
import sys
import pickle
import json
import itertools as it
from pprint import pprint
import collections as coll

import click

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from genome_mapping import utils
from genome_mapping import mappers
from genome_mapping import matchers
from genome_mapping import formatters
from genome_mapping.intervals import Tree
from genome_mapping.data import RESULT_TYPE


class ReadablePickleFile(click.File):
    name = 'readable-pickle'

    def __init__(self):
        super(ReadablePickleFile, self).__init__(mode='rb')

    def convert(self, value, param, ctx):
        fileobj = super(ReadablePickleFile, self).convert(value, param, ctx)
        return pickle.load(fileobj)


class WritablePickleFile(click.File):
    name = 'writeable-pickle'

    def __init__(self):
        super(WritablePickleFile, self).__init__(mode='wb')

    def convert(self, value, param, ctx):
        fileobj = super(WritablePickleFile, self).convert(value, param, ctx)
        return lambda d: pickle.dump(d, fileobj)


class KeyValue(click.ParamType):
    name = 'key-value'

    def convert(self, value, param, ctx):
        if value is None or value == {}:
            return value
        key, value = value.split('=')
        return {key: value}


@click.group()
def cli():
    pass


@cli.command('find')
@click.argument('genome', type=click.Path(readable=True))
@click.argument('targets', type=click.Path(exists=True, readable=True))
@click.argument('save', type=WritablePickleFile())
@click.option('--method', default='blat',
              type=click.Choice(mappers.known()))
@click.option('--organism', default='UNKNOWN')
def find(genome, targets, save, method='blat', organism='UNKNOWN'):
    """
    Search the genome for the given targets using the specified program.

    Parameters
    ----------

    genome :
        The path to the FASTA file of the genome to search in.

    targets :
        The path to the FASTA file of targets to search with.

    save :
        The path to the file to save pickled objects to.
    """
    mapper_class = mappers.fetch(method)
    mapper = mapper_class()
    save(list(mapper(genome, targets)))


@cli.group('hits')
def hits():
    """
    Group of commands dealing with the hits.
    """
    pass


@hits.command('select')
@click.argument('hits', type=ReadablePickleFile())
@click.argument('matcher', type=click.Choice(matchers.known()))
@click.argument('save', type=WritablePickleFile())
@click.option('--define', multiple=True, default={}, type=KeyValue())
def hits_select(hits, matcher, save, define={}):
    """
    Filter out the hits that only pass a match criteria. A match criteria is
    something like 'exact' or sequence identity.
    """
    definitions = {}
    for definition in define:
        definitions.update(definition)

    matcher_class = matchers.fetch(matcher)
    matcher = matcher_class(**definitions)
    save(list(matcher.filter_matches(hits)))


@hits.command('select-using-spec')
@click.argument('hits', type=ReadablePickleFile())
@click.argument('spec-file', type=click.File('rb'))
@click.argument('save', type=WritablePickleFile())
def hits_select_spec(hits, spec_file, save):
    """
    Select hits using the specifications in the given file. The file be a json
    object with a matcher entry that is the name of the matcher to use. It may
    also contain a JSON object of definitions to build the matcher with.
    """
    spec = json.load(spec_file)
    matcher = spec['matcher']
    if matcher not in matchers.known():
        raise ValueError("Unknown Matcher")
    definitions = spec.get('definitions', {})
    matcher_class = matchers.fetch(matcher)
    matcher = matcher_class(**definitions)
    save(list(matcher.filter_matches(hits)))


@hits.command('compare')
@click.argument('hits', type=ReadablePickleFile())
@click.argument('correct', type=click.Path(exists=True, readable=True))
@click.argument('save', type=WritablePickleFile())
def compare_matches(hits, correct, save):
    """
    Compare some hits to known examples to see how well they overlap.
    """
    tree = Tree(correct)
    save(tree.compare_to_known(hits))


@hits.command('best-within')
@click.argument('hits', type=ReadablePickleFile())
@click.argument('features', type=click.Path(exists=True, readable=True))
@click.argument('max-range', type=int)
@click.argument('save', type=WritablePickleFile())
def best_within(hits, features, max_range, save):
    """
    Find the best hits within some distance of the given features.
    """
    tree = Tree(features)
    save(tree.best_hits_within(hits, max_range))


@hits.command('merge')
@click.argument('hits', type=ReadablePickleFile(), nargs=-1)
@click.argument('save', type=WritablePickleFile())
def merge_hits(hits, save):
    """Merge several collections hits. This should not produce any duplicate
    hits.
    """
    save(sorted(set(it.chain.from_iterable(hits))))


@cli.group('comparisons')
def comparisons():
    """Group of commands dealing with maninpulating comparisons.
    """
    pass


@comparisons.command('select')
@click.argument('comparisons', type=ReadablePickleFile())
@click.argument('filter', type=str, nargs=-1)
@click.argument('save', type=WritablePickleFile())
def comparisons_select(comparisons, filter, save):
    """
    Filter comparisons to only those of the given type(s).
    """

    replace = {'is': '=='}
    processed = ' '.join(replace.get(w, w) for w in filter)
    ast = compile(processed, '<string>', mode='eval')

    def checker(obj):
        def hashable(name):
            return isinstance(getattr(obj, name), coll.Hashable)

        fields = utils.properities_of(obj.__class__)
        locals = {(f, getattr(obj, f)) for f in fields if hashable(f)}
        each = it.chain.from_iterable(t.split('_') for t in RESULT_TYPE)
        locals.update((t, t) for t in it.chain(RESULT_TYPE, each))
        return eval(ast, {}, dict(locals))

    save([comp for comp in comparisons if checker(comp)])


@comparisons.command('extract')
@click.argument('comparisons', type=ReadablePickleFile())
@click.argument('property')
@click.argument('save', type=WritablePickleFile())
@click.option('--skip-missing', is_flag=True, default=False)
def comparisons_extract(comparisons, property, save, skip_missing=False):
    """Extract parts of the given comparisons.
    """
    extracted = []
    for comparision in comparisons:
        entry = getattr(comparision, property)
        if skip_missing and not entry:
            continue
        extracted.append(entry)
    save(extracted)


@comparisons.command('summary')
@click.argument('comparisons', type=ReadablePickleFile())
@click.argument('save', type=click.File(mode='wb'))
def comparisons_summarize(comparisons, save):
    """
    Compute a summary of the number of each type of comparisons.
    """
    summary = dict(coll.Counter(c.type.pretty for c in comparisons))
    writer = csv.DictWriter(save, sorted(summary.keys()))
    writer.writeheader()
    writer.writerow(summary)


@comparisons.command('detect-variants')
@click.argument('comparisons', type=ReadablePickleFile())
@click.argument('save', type=WritablePickleFile())
def comparisons_detect_variants(comparisons, save):
    """Detect an splicing hits that are due to disagreement in locations
    between splicing variants. For example a location may have the splice
    variant noted, but the unspliced variant matches there as well.
    """
    # For all incorrect_exact matches
    # Find what the correct urs in the location is
    # See if there is a correct hit in the location
    # If in_hit has > 1 fragment and the correct hit has 1 => Spliced match transcript
    # If in_hit has 1 fragment and the correct has > 1 => transcript match spliced
    # Else leave alone
    pass


@comparisons.group('group')
def comparisons_group():
    pass


@comparisons_group.command('by-hit-urs')
@click.argument('comparisons', type=ReadablePickleFile())
@click.argument('save', type=WritablePickleFile())
def comparisons_group_hit_urs(comparisons, save):
    def key(comparision):
        if comparision.hit:
            return comparision.hit.urs
        return None
    ordered = sorted(comparisons, key=key)
    groups = it.groupby(ordered, key)
    save({urs: list(comps) for urs, comps in groups})


@comparisons_group.command('summarize-types')
@click.argument('grouped', type=ReadablePickleFile())
@click.argument('save', type=click.File(mode='wb'))
def comparisons_group_summary_type(grouped, save):
    header = ['urs'] + list(RESULT_TYPE)
    writer = csv.DictWriter(save, header)
    writer.writeheader()
    for urs, hits in grouped.iteritems():
        entry = coll.defaultdict(int)
        entry['urs'] = urs
        entry.update(coll.Counter(h.type.pretty for h in hits))
        writer.writerow(entry)


@cli.command('as')
@click.argument('format', type=click.Choice(formatters.known()))
@click.argument('data', type=ReadablePickleFile())
@click.argument('save', type=click.File(mode='wb'))
def format(format, data, save):
    """
    Format the data into some format.

    Paramaters
    ----------
    format : str
        The name of the format to use.
    data : path
        The path to the data to read, '-' means stdin.
    save : path
         The path of where to save data, '-' means stdout.
    """
    formatters.format(data, format, save)


@cli.command('pp')
@click.argument('data', type=ReadablePickleFile())
def display(data):
    """
    Pretty print data to stdout.
    """
    pprint(data)


@cli.command('head')
@click.argument('data', type=ReadablePickleFile())
@click.argument('save', type=WritablePickleFile())
@click.option('--count', type=int, default=10)
def head(data, save, count):
    return slice_of(data, 0, count, save)


# @cli.command('slice')
# @click.argument('data', type=ReadablePickleFile())
# @click.argument('start', type=int, default=0)
# @click.argument('stop', type=int, default=10)
# @click.argument('save', type=WritablePickleFile())
def slice_of(data, start, stop, save):
    if isinstance(data, (tuple, list)):
        return save(data[start:stop])
    if isinstance(data, dict):
        keys = sorted(data.iterkeys())[start:stop]
        return save({k: data[k] for k in keys})
    raise ValueError("Cannot slice given type")


if __name__ == '__main__':
    cli()
