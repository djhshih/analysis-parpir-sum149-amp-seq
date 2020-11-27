library(io)
library(dplyr)

pheno <- qread("../annot/pheno.tsv");

biosample.d <- transmute(pheno,
	sample_name = id,
	sample_title = description,
	organism = "Homo sapiens",
	isolate = paste0(id, "_b", batch),
	age = ifelse(cell_line == "SUM149PT", 40, 36),
	biomaterial_provider = ifelse(cell_line == "SUM149PT", "BioIVT", "ATCC"),
	sex = "female",
	tissue = "mammary gland",
	cell_line = cell_line
);

qwrite(biosample.d, "biosample-attributes.tsv")

sra.d <- transmute(pheno,
	sample_name = id,
	library_ID = id,
	title = paste0("Amplicon-Seq of BRCA1 exon 10 locus: SUM149PT"),
	library_strategy = "AMPLICON",
	library_source = "GENOMIC",
	library_selection = "PCR",
	library_layout = "paired",
	platform = "ILLUMINA",
	instrument_model = "Illumina MiSeq",
	design_description = "Amplicon-EZ",
	filetype = "fastq",
);

# pair fastq files together

files <- list_files("../fastq", pattern=".+fastq.gz");

files.d <- data.frame(
	sample_id = sub("_.+", "", sub("[^-]+-", "", files)),
	filename = files,
	read = sub(".*(R\\d).*", "\\1", files)
);

r1.d <- filter(files.d, read == "R1");
r2.d <- filter(files.d, read == "R2");

paired.d <- left_join(
	select(r1.d, -read),
	select(r2.d, -read),
	suffix = c("", "2"),
	by = "sample_id"
);

sra.d <- left_join(sra.d, paired.d, by = c(sample_name="sample_id"));

qwrite(sra.d, "sra-metadata.tsv");

