#!/usr/bin/env python

import csv
import click
import gffutils


class GffInput(click.File):
    def __init__(self):
        super(GffInput, self).__init__(mode='rb')

    def convert(self, value, param, ctx):
        with super(GffInput, self).convert(value, param, ctx) as fileobj:
            return gffutils.create_db(fileobj.read(), ':memory:',
                                      from_string=True)


@click.command()
@click.argument('gff', type=GffInput())
@click.argument('bed', type=click.File(mode='wb'))
def main(gff, bed):
    writer = csv.writer(bed, delimiter='\t', quoting=csv.QUOTE_NONE)
    for transcript in gff.features_of_type('transcript'):
        start = transcript.start
        exons = list(gff.children(transcript.id))
        if len(exons) == 0:
            raise ValueError("No exons found")

        exon_sizes = ','.join(str(e.end - e.start + 1) for e in exons)
        exon_starts = ','.join(str(start - e.start) for e in exons)

        writer.writerow([
            transcript.seqid,
            start - 1,
            transcript.end,
            transcript.attributes['Name'][0],
            0,
            transcript.strand,
            start,
            transcript.end,
            '63,125,151',
            len(exons),
            exon_sizes,
            exon_starts,
        ])


if __name__ == '__main__':
    main()
