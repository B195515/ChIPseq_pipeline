# Generate QC information of short reads
library(ShortRead)
fls <- list.files(pattern="*.fastq*", full=TRUE)
names(fls) <- sub(".fastq", "", basename(fls))
qas <- lapply(seq_along(fls), 
	function(i, fls) qa(
		readFastq(fls[i]),
		names(fls)[i],
		fls))
qa <- do.call(rbind, qas)
rpt <- report(qa, dest='QA_report.html')

# Build annotation files
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
si <- seqinfo(genome)
si <- si[paste0('chr', c(1:19,'X','Y'))]

# Get mart from ENSEMBL for annotation
library(biomaRt)
mart <- useMart(
	biomart="ENSEMBL_MART_ENSEMBL",
	dataset="mmusculus_gene_ensembl",
	host="www.ensembl.org")
fm <- Gviz:::.getBMFeatureMap()


# Get snapshot of Chr6 from 
# pos 122530000 to pos 122900000
# encoding an ES-specific Nanog gene
# Isolate gene models for this interval
bm <- BiomartGeneRegionTrack(
	chromosome='chr6',
	genome="mm10",
	start=122530000,
	end=122900000,
	biomart=mart,
	filter=list("with_refseq_mrna"=TRUE),
	size=4,
	name="RefSeq",
	utr5="red3", utr3="red3",
	protein_coding="black",
	col.line=NULL,
	cex=7,
	collapseTranscripts="longest",
	featureMap=fm)

# Build a basic annotation resource
listAttributes(mart)[1:3,]
ds <- useDataset('mmusculus_gene_ensembl', mart=mart)
chroms <- 6
egs <- getBM(
	attributes=c('ensembl_gene_id', 
		'external_gene_name', 'chromosome_name',
		'start_position', 'end_position', 'strand'),
	filters='chromosome_name',
	values=chroms,
	mart=ds)

#------------(Shell: MACS14 peak calling - filtered Chr6 only)

# Start actual analysis

library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(chipseq)
load("chip_seq.Rdata")

# load BED files before MACS14
input <- import.bed(file.path(getwd(), 'ES_input_filtered_ucsc_chr6.bed'))
rep1 <- import.bed(file.path(getwd(), 'H3K27ac_rep1_filtered_ucsc_chr6.bed'))
rep2 <- import.bed(file.path(getwd(), 'H3K27ac_rep2_filtered_ucsc_chr6.bed'))
# check similar lengths
length(input)
length(rep1)
length(rep2)

# Extend reads by median fragment length
# Want peaks at fragment centres
prepareChIPseq <- function(reads){
	frag.len <- median(estimate.mean.fraglen(reads))
	cat(paste0('Median fragment size for this library is ',
		round(frag.len)))
	reads.extended <- resize(reads, width=frag.len)
	return(trim(reads.extended))
}
# Run the function for each sample
# These are GRanges object with 6 ranges and 2 metadata columns
input <- prepareChIPseq(input)
rep1 <- prepareChIPseq(rep1)
rep2 <- prepareChIPseq(rep2)
# Check samples
head(rep1)

# Count reads mapped to predefine genomic intervals
# build bins from si annotation object
# also GRanges object with 6 ranges
binsize <- 200
bins <- tileGenome(
	si['chr6'],
	tilewidth=binsize,
	cut.last.tile.in.chrom=TRUE)
head(bins)

# Count elements that fall into each bin
BinChIPseq <- function(reads,bins){
	mcols(bins)$score <- countOverlaps(bins,reads)
	return(bins)
}
input.200bins <- BinChIPseq(input,bins)
rep1.200bins <- BinChIPseq(rep1,bins)
rep2.200bins <- BinChIPseq(rep2,bins)
# Check output
head(rep1.200bins)

#*Plot coverage for 1000 bins from bin 200k
png("plot_coverage_bin200k_chr6.png")
plot(200000:201000,
	rep1.200bins$score[200000:201000],
	xlab='chr6',
	ylab='counts per bin',
	main='Coverage for 1000 bins from bin 200k on Chromosome 6')
dev.off()

# Export chipseq data as bedGraph
# To visualized in UCSC viewer
export(input.200bins,
	con='input_chr6.bedGraph',
	format='bedGraph')
export(rep1.200bins,
	con='H3K27ac_rep1_chr6.bedGraph',
	format='bedGraph')
export(rep2.200bins,
	con='H3K27ac_rep2_chr6.bedGraph',
	format='bedGraph')

# Visualize tracks using GViz library
# Using bm annotation object
library(Gviz)
# create a genomic axis track
AT <- GenomeAxisTrack()
# use plot tracks to plot a genomic region
plotTracks(c(bm,AT),
	from=122530000,
	to=122900000,
	transcriptAnnotation="symbol",
	window="auto",
	cex.title=1,
	fontsize=10)
# build extra tracks for visualization
input.track <- DataTrack(
	input.200bins,
	strand="*",
	genome="mm10",
	col.histogram="gray",
	fill.histogram="black",
	name="Input",
	col.axis="black",
	cex.axis=0.4,
	ylim=c(0,150))
rep1.track <- DataTrack(
	rep1.200bins,
	strand="*",
	genome="mm10",
	col.histogram="steelblue",
	fill.histogram="black",
	name="Rep.1",
	col.axis="steelblue",
	cex.axis=0.4,
	ylim=c(0,150))
rep2.track <- DataTrack(
	rep2.200bins,
	strand="*",
	genome="mm10",
	col.histogram="steelblue",
	fill.histogram="black",
	name="Rep.2",
	col.axis="steelblue",
	cex.axis=0.4,
	ylim=c(0,150))

# Now plot the results
#TODO: png resolution?
png("plottracks_input.rep1.rep2.bm.AT.png",
	width=1200, height=1200)
plotTracks(
	c(input.track, rep1.track, rep2.track, bm, AT),
	from=122530000,
	to=122900000,
	transcriptAnnotation="symbol",
	window="auto",
	type="histogram",
	cex.title=0.7,
	fontsize=10)
dev.off()


#-------Analysis of peaks from MACS software

# Load MACS peaks into R
peaks.rep1 <- import.bed(file.path(getwd(), 'Rep1_peaks_ucsc_chr6.bed'))
peaks.rep2 <- import.bed(file.path(getwd(), 'Rep2_peaks_ucsc_chr6.bed'))

#Create annotation tracks from Peaks
peaks1.track <- AnnotationTrack(
	peaks.rep1,
	genome="mm10",
	name="Peaks Rep. 1",
	chromosome="chr6",
	shape="box",
	fill="blue3",
	size=2)
peaks2.track <- AnnotationTrack(
	peaks.rep2,
	genome="mm10",
	name="Peaks Rep. 2",
	chromosome="chr6",
	shape="box",
	fill="blue3",
	size=2)
# Plot the peak tracks - focusing on a region near the Nanog gene
png("plottrack_input.rep1.rep2.peaks1.peaks2.bm.AT.png",
	width=1000, height=1000)
plotTracks(
	c(input.track, rep1.track, peaks1.track, 
		rep2.track, peaks2.track, bm, AT),
	from=122700000,
	to=122730000,
	transcriptAnnotation="symbol",
	window="auto",
	type="histogram",
	cex.title=0.7,
	fontisze=10)
dev.off()

# Venn Diagrams to show overlap between Rep1 and Rep2
ovlp <- findOverlaps(peaks.rep1, peaks.rep2)
ovlp
# number of overlaps that one peak in sample 1 has to 
# multiple peaks in sample 2
ov <- min(
	length(unique(queryHits(ovlp))),
	length(unique(subjectHits(ovlp))))
ov
# Create the Venn Diagram
library(VennDiagram)
png("venn_peaksrep1.peaksrep2.ov.png",
	width=1000, height=1000)
draw.pairwise.venn(
	area1=length(peaks.rep1),
	area2=length(peaks.rep2),
	cross.area=ov,
	category=c("rep1", "rep2"),
	fill=c("steelblue", "blue3"),
	cat.cex=0.7)
dev.off()


# Peaks identified in both reps (enriched regions)
# Id enriched regions
enriched.regions <- Reduce(
	subsetByOverlaps, 
	list(peaks.rep1, peaks.rep2))
# Build annotation track for the enriched regions
enr.reg.track <- AnnotationTrack(
	enriched.regions,
	genome="mm10",
	name="Enriched regions",
	chromosome="chr6",
	shape="box",
	fill="green3",
	size=2)
# Plot the tracks again with the extra info
png("plottracks_input.rep1.rep2.peaks1.peaks2.enrreg.bm.AT.png",
	width=1000, height=1000)
plotTracks(
	c(input.track, rep1.track, peaks1.track,
		rep2.track, peaks2.track, enr.reg.track,
		bm, AT),
	from=122700000,
	to=122750000,
	transcriptAnnotation="symbol",
	window="auto",
	type="histogram",
	cex.title=0.5,
	fontsize=16)
dev.off()


#------Isolation of promoters overlapping H3K27ac peaks

# Identify the TSS, consider gene orientation
egs$TSS <- ifelse(
	egs$strand=="1",
	egs$start_position,
	egs$end_position)
head(egs)

# if promoters at TSS region +/-200bp
promoter_regions <- GRanges(
	seqnames=Rle(paste0("chr", egs$chromosome_name)),
	ranges=IRanges(start=egs$TSS-200,
		end=egs$TSS+200),
	strand=Rle(rep("*", nrow(egs))),
	gene=egs$external_gene_id)
promoter_regions

#------Overlap promoters with Enriched Regions
ovlp2b <- findOverlaps(promoter_regions, enriched.regions)
cat(sprintf(
	"%d of %d enriched regions overlap a promoter",
	length(unique(subjectHits(ovlp2b))), length(enriched.regions)
	))
promoter_total_length <- sum(width(reduce(promoter_regions)))
promoter_total_length

# Fraction of the chromosome
promoter_frac_of_chr6 <- promoter_total_length / seqlengths(si)["chr6"]
promoter_frac_of_chr6

# Confirm enrichment
binom.enr <- binom.test(
	length(unique(subjectHits(ovlp2b))),
	length(enriched.regions),
	promoter_frac_of_chr6
	)
# The estimated probability of success
binom.enr$estimate

# Which prom overlap with H327ac peaks?
pos.TSS <- egs[unique(queryHits(findOverlaps(
	promoter_regions, enriched.regions)
)),]
pos.TSS[1:3,]


#-----H3K27ac distribution around a subset of gene proms
# Extend prom regions to +/-1000bp around the TSS
# Split regions into 20 bins
# Order regions by gene strand
# Tile regions with these bins, consider gene orientation
tiles <- sapply(1:nrow(pos.TSS), function(i)
	if(pos.TSS$strand[i]=="1")
		pos.TSS$TSS[i]+seq(-1000,900,length.out=20)
	else
		pos.TSS$TSS[i]+seq(900,-1000,length.out=20)
		)

tiles <- GRanges(
	tilename=paste(
		rep(pos.TSS$ensembl_gene_id, each=20),
		1:20, sep="_"),
	seqnames=Rle(rep(
		paste0("chr", pos.TSS$chromosome_name),
		each=20)),
	ranges=IRanges(start=as.vector(tiles), width=100),
	strand=Rle(rep("*", length(as.vector(tiles)))),
	seqinfo=si
	)

tiles

# Count reads that map to the tile
# Generate a matrix where
# each row reps an enriched promoter
# each col reps promoter tiles sequence
H3K27ac.p <- countOverlaps(tiles,rep1)+countOverlaps(tiles,rep2)
H3K27ac.p.matrix <- matrix(
	H3K27ac.p, 
	nrow=nrow(pos.TSS),
	ncol=20,
	byrow=TRUE)


# Plot the results
# layout and par set the layout and margins
# Other commands: Box, Image, Plot
png("H3K27ac_distribution.png",
	width=1000, height=1000)

colors <- colorRampPalette(c('white','red','gray','black'))(100)

layout(
	mat=matrix(c(1,2,0,3),2,2),
	widths=c(2,2,2),
	heights=c(0.5,5,0.5,5),
	TRUE)
par(mar=c(2,2,1.5,1))

# plot the legend
image(
	seq(0, max(H3K27ac.p.matrix), length.out=100),
	1,
	matrix(seq(0, max(H3K27ac.p.matrix), length.out=100),
		100,1),
	col=colors,
	xlab="Distance from TSS",
	ylab="",
	main="Number of reads",
	yaxt="n",
	lwd=3,
	axes=TRUE)
box(col="black", lwd=2)

# left side
image(
	x=seq(-1000,1000, length.out=20),
	y=1:nrow(H3K27ac.p.matrix),
	z=t(H3K27ac.p.matrix[order(rowSums(H3K27ac.p.matrix)),]),
	col=colors,
	xlab="Distance from TSS (bp)",
	ylab="Promoters",
	lwd=2)
box(col="black", lwd=2)
abline(v=0, lwd=1, col="gray")

# right side - scatter plot
plot(
	x=seq(-1000,1000,length.out=20),
	y=colMeans(H3K27ac.p.matrix),
	ty="b", pch=19,
	col="red4", lwd=2,
	ylab="Mean tag count",
	xlab="Distance from TSS (bp)")
#background grids
abline(
	h=seq(1,100,by=5),
	v=seq(-1000,1000,length.out=20),
	lwd=0.25, col="gray")
# borders
box(col="black", lwd=2)
dev.off()