gm=bin/gm.py
all_md5=data/md5.tsv

all : pombe zebrafish worm fly

###############################################################################
# pombe
###############################################################################

pombe : data/pombe/inferred.json data/pombe/summary.csv data/pombe/unknown.fasta
	@echo "POMBE: Found `jq 'length' data/pombe/inferred.json` coordinates of `grep -c '^>' data/pombe/unknown.fasta`"

data/pombe/genome.fasta : bin/fetch
	$^ "ftp://ftp.ensemblgenomes.org/pub/release-35/fungi/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz" > $@

data/pombe/known.gff3 : bin/fetch
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates/Schizosaccharomyces_pombe.ASM294v2.gff3.gz" > $@

data/pombe/raw.bed : bin/fetch-bed
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates/Schizosaccharomyces_pombe.ASM294v2.bed.gz" > $@


###############################################################################
# C. elegans
###############################################################################

worm : data/worm/inferred.json data/worm/summary.csv data/worm/unknown.fasta
	@echo "WORM: Found `jq 'length' data/worm/inferred.json` coordinates of `grep -c '^>' data/worm/unknown.fasta`"

data/worm/genome.fasta : bin/fetch
	$^ "ftp://ftp.wormbase.org/pub/wormbase/releases/WS251/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS251.genomic.fa.gz" > $@

data/worm/known.gff3 : bin/fetch
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates/Caenorhabditis_elegans.WBcel235.gff3.gz" > $@

data/worm/raw.bed : bin/fetch-bed
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates/Caenorhabditis_elegans.WBcel235.bed.gz" > $@

###############################################################################
# Drosophila melanogaster (6,515)
###############################################################################

fly: data/fly/inferred.json data/fly/summary.csv data/fly/unknown.fasta
	@echo "FLY: Found `jq 'length' data/fly/inferred.json` coordinates of `grep -c '^>' data/fly/unknown.fasta`"

data/fly/genome.fasta : bin/fetch
	$^ "ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.15_FB2017_02/fasta/dmel-all-chromosome-r6.15.fasta.gz" > $@

data/fly/known.gff3 : bin/fetch
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates//Drosophila_melanogaster.BDGP6.gff3.gz" > $@

data/fly/raw.bed : bin/fetch-bed
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates//Drosophila_melanogaster.BDGP6.bed.gz" > $@


###############################################################################
# Drosophila melanogaster (6,515)
###############################################################################

zebrafish: data/zebrafish/inferred.json data/zebrafish/summary.csv data/zebrafish/unknown.fasta
	@echo "ZEBRAFISH: Found `jq 'length' data/zebrafish/inferred.json` coordinates of `grep -c '^>' data/zebrafish/unknown.fasta`"

data/zebrafish/genome.fasta : bin/fetch
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.1.fa.gz" > $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.2.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.3.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.4.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.5.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.6.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.7.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.8.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.9.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.10.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.11.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.12.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.13.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.14.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.15.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.16.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.17.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.18.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.19.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.20.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.21.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.22.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.23.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.24.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.25.fa.gz" >> $@
	$^ "ftp://ftp.ensembl.org/pub/release-88/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.chromosome.MT.fa.gz" >> $@

data/zebrafish/known.gff3 : bin/fetch
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates//Danio_rerio.GRCz10.gff3.gz" > $@

data/zebrafish/raw.bed : bin/fetch-bed
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/6.0/genome_coordinates//Danio_rerio.GRCz10.bed.gz" > $@


###############################################################################
# Generic Definitions
###############################################################################

$(all_md5) : bin/fetch
	$^ "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral//releases/6.0/md5/md5.tsv.gz" > $@
	xsv index $@

%.2bit : %.fasta
	faToTwoBit $^ $@

%.11ooc : %.fasta
	blat -makeOoc=11.ooc $^

data/%/known-md5.txt : bin/known-md5s $(all_md5) data/%/known.gff3 data/%/bad-ids
	$^ > $@

data/%/targets.bed : bin/targets-bed data/%/raw.bed data/%/known.gff3 data/%/bad-ids
	$^ > $@

data/%/targets.fasta : bin/targets data/%/genome.fasta data/%/targets.bed data/%/md5sums.txt data/%/known-md5.txt
	$^ > $@

data/%/targets-hits.pickle : $(gm) data/%/genome.2bit data/%/targets.fasta
	date
	$(gm) find $(word 2,$^) $(word 3,$^) $@
	date

data/%/unknown-hits.pickle : $(gm) data/%/genome.2bit data/%/unknown.fasta
	$(gm) find $(word 2,$^) $(word 3,$^) $@

data/%/targets-compared.pickle : $(gm) data/%/targets-selected.pickle data/%/known.gff3
	$(gm) hits compare $(word 2,$^) $(word 3,$^) $@

data/%/compared.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) as gff3 $(word 2,$^) $@

data/%/missing.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons select $(word 2,$^) type.pretty is missing - | $(gm).py as gff3 - $@

data/%/novel.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons select $(word 2,$^) type.pretty is novel - | $(gm).py as gff3 - $@

data/%/incorrect.gff3 : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons select $(word 2,$^) type.pretty is incorrect - | $(gm).py as gff3 - $@

data/%/summary.csv : $(gm) data/%/targets-compared.pickle
	$(gm) comparisons summary $(word 2,$^) - | tee $@ | xsv table

data/%/targets-selected.pickle : $(gm) data/%/targets-hits.pickle data/%/select-spec.json
	$(gm) hits select-using-spec $(word 2,$^) $(word 3,$^) $@

data/%/unknown-selected.pickle : $(gm) data/%/unknown-hits.pickle data/%/select-spec.json
	$(gm) hits select-using-spec $(word 2,$^) $(word 3,$^) $@

data/%/inferred.json : $(gm) data/%/unknown-selected.pickle
	$(gm) as insertable $(word 2,$^) $@

.PHONY : all
