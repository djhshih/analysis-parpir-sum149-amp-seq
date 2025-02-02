library(ggplot2)
library(io)

# Calculate the probability of having at least one revertant among n genomes,
# given p, the probability of having a revertant mutation in one genome
revertant_prob <- function(n, p) {
	1 - (1 - p)^n
}

# Zamborszky 2016 https://doi.org/10.1038/onc.2016.243
# Genomes of BRCA1−/− mock-treated samples contained 8.0±1.0 short insertions and 12.7±1.2 deletions
# BRCA2−/− mock-treated samples contained insertions (10.3±3.1, P=0.045) and deletions (40.3±2.1, P<0.001).
# There was no significant difference between the numbers of spontaneous short indels in the WT sample
# and the BRCA1+/− and BRCA2+/− heterozygotes

# BRCA1 has 1884 aa = 5652 bp CDS 
# BRCA2 has 3418 aa = 10254 bp CDS

# GRCh38.p13 has 3,110,748,599 non-N bases

# second indel can be before or after the original frameshift mutation
# the limiting factor is where the first indel introduces a pre-mature stop
# codon due to the frameshift, as well as where the second indel would
# introduce a pre-mature stop codong due to frameshift


# number of cells
n <- 1e6;

rescuable <- 120;
rescuable.lower <- 60;
rescuable.upper <- 600;

# BRCA1

target <- 5652;
genome <- 3110748599;

# restore a 1-bp deletion

# expected number of 1-bp insertion per genome
d <- 5;

# REMARK: To consider variations in d, we can sample d
#         from a distribution such as the Poisson distribution.

ns <- 10^seq(1, 8, 0.05);

# probability of having a revertant mutation in one genome
p <- dbinom(1, d, rescuable / genome);
probs <- vapply(ns, revertant_prob, p = p, 0);

revertant_prob(1e6, p)

p.lower <- dbinom(1, d, rescuable.lower / genome);
probs.lower <- vapply(ns, revertant_prob, p = p.lower, 0);

p.upper <- dbinom(1, d, rescuable.upper / genome);
probs.upper <- vapply(ns, revertant_prob, p = p.upper, 0);

d.brca1 <- data.frame(
	cells = ns,
	prob = probs,
	prob_lower = probs.lower,
	prob_upper = probs.upper
);

qdraw(
	ggplot(d.brca1, aes(x = cells, y = prob, ymin = prob_lower, ymax = prob_upper)) + 
		theme_bw() +
		geom_ribbon(alpha=0.2, fill="firebrick") +
		geom_line() + 
		scale_x_log10() + annotation_logticks(sides="b") +
		xlab("number of cancer cells") + ylab("revertant probability") +
		ggtitle("BRCA1")
	,
	file = "revertant-prob_brca1.pdf"
)


# BRCA2

genome <- 3110748599;

# expected number of 1-bp insertion per genome
d <- 6;

revertant_prob(1e6, dbinom(1, d, rescuable / genome))
revertant_prob(1e6, dbinom(1, 1, rescuable / genome))

# probability of having a revertant mutation in one genome
p <- dbinom(1, d, rescuable / genome);
probs <- vapply(ns, revertant_prob, p = p, 0);

p.lower <- dbinom(1, d, rescuable.lower / genome);
probs.lower <- vapply(ns, revertant_prob, p = p.lower, 0);

p.upper <- dbinom(1, d, rescuable.upper / genome);
probs.upper <- vapply(ns, revertant_prob, p = p.upper, 0);

d.brca2 <- data.frame(
	cells = ns,
	prob = probs,
	prob_lower = probs.lower,
	prob_upper = probs.upper
);

qdraw(
	ggplot(d.brca2, aes(x = cells, y = prob, ymin = prob_lower, ymax = prob_upper)) + 
		theme_bw() +
		geom_ribbon(alpha=0.2, fill="firebrick") +
		geom_line() + 
		scale_x_log10() + annotation_logticks(sides="b") +
		xlab("number of cancer cells") + ylab("revertant probability") +
		ggtitle("BRCA2")
	,
	file = "revertant-prob_brca2.pdf"
)

