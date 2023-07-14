library(snpR)
d <- readRDS("Data/monarch_nomaf.RDS")
d <- d[pop = c("NAM", "ROT", "GUA", "HAW", "NSW", "QLD", "ROT", "SAI")]
d <- filter_snps(d, mgc = 1, hwe = 1e-6, hwe_facets = "pop")

set.seed(1334512)
sub <- sort(sample(nrow(d), 5e4, replace = FALSE))

dsub <- d[sub,]
nm <- sample.meta(dsub)
nm$pop[nm$pop == "NAM"] <- ifelse(grepl("Mexico", nm$samp[which(nm$pop == "NAM")]), "ENA", "WNA")
nm$pop[nm$samp == "NAM_M9.10"] <- "ENA"
colnames(dsub) <- paste0(nm$pop, 1:nrow(nm))
sample.meta(dsub) <- nm
format_snps(dsub, "vcf", outfile = "data/monarchs.vcf", chr = "group")
