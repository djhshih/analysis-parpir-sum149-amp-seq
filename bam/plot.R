library(io)
library(ggsci)
library(ggplot2)
library(dplyr)
library(tidyr)
library(binom)  # binom.confint

logistic <- function(x) {
	1 / (1 + exp(-x))
}

logit <- function(x) {
	log(x) - log(1 - x)
}

my_theme <- function() {
	theme_bw() +
	theme(
		panel.grid = element_blank()
	)
}

binom_confint <- function(x, n) {
	binom.confint(x, n, method="exact")
}

####

out.fn <- filename("brca-amp-seq");
pdf.fn <- insert(out.fn, ext="pdf");

d <- qread("stats.tsv");

d$total <- d$wildtype + d$parental + d$revertant;

####

# estimate adapter cross contamination in batch A

da <- d[d$batch == "A", ];
db <- d[d$batch == "B", ];
mcf10a <- which(da$cell_line == "MCF10A");

# independent

# MCF10A should not have parental and revertant reads
# SUM149 should not have wildtype types
f.ind <- c(
	# MCF10A only
	(da$parental[mcf10a] + da$revertant[mcf10a]) / da$total[mcf10a],
	# SUM149 cell lines
	da$wildtype[-mcf10a] / da$total[-mcf10a]
);

f.mean <- mean(f.ind);
f.mean.ci.width <- qnorm(1 - 0.025) * sd(f.ind) / sqrt(length(f.ind));
c(f.mean, f.mean - f.mean.ci.width, f.mean + f.mean.ci.width)
# 0.40% (0.16% - 0.64%)

# independent, in log scale

lf.ind <- logit(f.ind);
lf.mean <- mean(lf.ind);
lf.mean.ci.width <- qnorm(1 - 0.025) * sd(lf.ind) / sqrt(length(lf.ind));
logistic( c(lf.mean, lf.mean - lf.mean.ci.width, lf.mean + lf.mean.ci.width) )
# 0.33% (0.14% - 0.75%)

# pooled

x <- (da$parental[mcf10a] + da$revertant[mcf10a] + sum(da$wildtype[-mcf10a]));
n <- (da$total[mcf10a] + sum(da$total[-mcf10a]));

# get binomial confidence interval
binom.test(x, n)
# 0.40% (0.38% - 0.42%)

####

# The three wildtype reads in DS02-P-3 look real.
# These could have come from low-level cross-contamination of the SUM149-P cell line
# DS01-P-2 and DS02-P-3 were processed independently from all other samples
# The possiblity that the parental cell was cross-contaminated by the
# resistant clones cannot be completely excluded, either.

# In batch A, use the count from the highest level of cross-contamination to 
# adjust for cross-contamination.
# This will undercount the reads in the source of the contamionation, but
# since the cross-contamination rate is low, this underestimation is
# negligible.
# Ideally, we would account for sample-specific read depth.

da$wildtype <- pmax(0, da$wildtype - max(da$wildtype[-mcf10a]));
da$parental <- pmax(0, da$parental - max(da$parental[mcf10a]));
da$revertant <- pmax(0, da$revertant - max(da$revertant[mcf10a]));

d.c <- rbind(da, db);

# $ faCount brca1-variants.fasta
# seq    len     A       C       G       T       N       cpg
# BRCA1_amplicon1 357     129     62      77      89      0       2
# BRCA1_amplicon1-parental        356     129     62      77      88      0	2
# BRCA1_amplicon1-revertant       321     112     56      69      84      0	2

# parental sequence is longer than revertant by 10.9%
356 / 321
# which likely causes a read depth bias factor of 1.37 - 1.47
bias.ind <- c(51597 / 37638, 61724 / 41754);
bias.mean <- mean(bias.ind);
bias.overall <- (51597 + 61724) / (37639 + 41754);

# pool samples together
d.s <- select(d.c, cell_line, wildtype, parental, revertant) %>% group_by(cell_line) %>%
	summarize(wildtype = sum(wildtype), parental = sum(parental), revertant = sum(revertant)) %>%
	ungroup();

d.s <- mutate(d.s, revertant = revertant / bias.overall);


parental.counts <- as.numeric(d.s[d.s$cell_line == "SUM149-P", c("parental", "revertant_adj")]);
revertant.cell.freq <- 2 * parental.counts[2] / sum(parental.counts);
# 2.43e-05

# expect this many cells to have one revertant cell
1 / revertant.cell.freq
# 1 in 41131


# normalize across rows
d.sn <- mutate(d.s, total = wildtype + parental + revertant) %>%
	group_by(cell_line) %>%
	summarize(
		wildtype = wildtype / total,
		parental = parental / total,
		revertant = revertant / total
	) %>% ungroup();

d.m <- pivot_longer(d.sn, c("wildtype", "parental", "revertant"), names_to = "allele")
d.m$allele <- factor(d.m$allele, c("wildtype", "parental", "revertant"));

qdraw(
	ggplot(d.m, aes(x=cell_line, y=value, fill=allele)) +
		geom_col(position="stack") + my_theme() +
		scale_fill_npg() +
		#coord_flip() + scale_x_discrete(limits=rev(levels(d.m$cell_line))) +
		theme(axis.text.x = element_text(angle=45, hjust=1)) +
		xlab("") + ylab("allelic frequency")
	,
	width = 3, height = 4,
	file = insert(pdf.fn, "vaf")
)

d.sn.revertant <- mutate(d.s, total = wildtype + parental + revertant) %>%
	group_by(cell_line) %>%
	summarize(
		y = revertant / total,
		ymin = pmax(0, binom_confint(revertant, total)$lower),
		ymax = pmin(1, binom_confint(revertant, total)$upper)
	) %>% ungroup();

d.sn.revertant[1, "ymax"] <- 0;

# coord_flip is not compatible with annotation_logticks

qdraw(
	ggplot(d.sn.revertant, aes(x=cell_line, y=y, ymin=ymin, ymax=ymax)) +
		geom_errorbar(width=0.1, colour="grey60") + my_theme() +
		geom_point() + 
		scale_y_log10(limits=c(1e-6, 0.75), breaks=c(0, 1e-2, 1e-5, 0.5, 1)) +
		annotation_logticks(sides="l", colour="grey60") +
		#coord_flip() + scale_x_discrete(limits=rev(levels(d.sn$cell_line))) +
		theme(axis.text.x = element_text(angle=45, hjust=1)) +
		xlab("") + ylab("frequency of revertant allele")
	,
	width = 2, height = 4,
	file = insert(pdf.fn, "revertant-vaf")
)

