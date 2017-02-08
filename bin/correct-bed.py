#!/usr/bin/env python

import csv

import click


@click.command()
@click.argument('bed', type=click.File(mode='rb'))
@click.argument('save', type=click.File(mode='wb'))
def main(bed, save):
    writer = csv.writer(save, delimiter='\t', quoting=csv.QUOTE_NONE)
    for row in csv.reader(bed, delimiter='\t'):
        row[1] = int(row[1]) - 1
        row[6] = int(row[6]) - 1
        row[10] = ','.join(str(int(size) + 1) for size in row[10].split(','))
        writer.writerow(row)


if __name__ == '__main__':
    main()
