coverage <- read.table("~/Bioinformática/TreinamentoAnaliseDNA/HCMV/aligned/coverage.txt", header=FALSE)
colnames(coverage) <- c("chrom","position","depth")

ggplot(coverage, aes(position, depth)) +
  geom_line() +
  theme_bw() +
  labs(
    x="Genomic position",
    y="Read depth",
    title="Genome coverage depth"
  )
