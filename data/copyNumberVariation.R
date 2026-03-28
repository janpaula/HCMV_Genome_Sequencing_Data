library(GenomicAlignments)
library(GenomicRanges)

# Importa o arquivo BAM
bam <- BamFile("~/Bioinformática/TreinamentoAnaliseDNA/HCMV/aligned/final_sorted.bam")

# Escaneia o arquivo BAM
scanBamHeader(bam) # Verifica numero do gene

seqlevels(genes) # Verifica numero do gene

# Altera para apenas um numero
seqlevels(genes) <- "NC_006273.2"
seqnames(genes) <- "NC_006273.2"

# Conta reads por gene
gene_counts <- summarizeOverlaps(
  features = genes,
  reads = bam,
  mode = "Union",
  singleEnd = TRUE,
  ignore.strand = TRUE
)

# Extrai contagens
counts <- assay(gene_counts)
length(counts)

# Verifica os nomes das colunas para fazer corretamente o cnv
colnames(mcols(genes))

# omite NAs
gene_name <- ifelse(
  is.na(genes$gene),
  genes$locus_tag,
  genes$gene
)

# CNV
cnv_df <- data.frame(
  gene = gene_name,
  start = start(genes),
  end = end(genes),
  counts = counts
)

# Renomeia "final_sorted.bam" por "counts
colnames(cnv_df)[4] <- "counts"

#Verifica se está ok
head(cnv_df)

# PASSO 2
# normalizar para identificar CNV
median_cov <- median(cnv_df$counts[cnv_df$counts > 0])

# Filtrar genes sem cobertura
cnv_df_filtered <- cnv_df[cnv_df$counts > 0, ]


# Verificar
nrow(cnv_df_filtered)

# Calculo da mediana
median_cov <- median(cnv_df_filtered$counts)

# Recalcular CNV
cnv_df_filtered$cnv_ratio <- cnv_df_filtered$counts / median_cov

# Verifica
summary(cnv_df_filtered$cnv_ratio)


# Grafico filtrado
plot(
  cnv_df_filtered$start,
  cnv_df_filtered$cnv_ratio,
  pch = 16,
  xlab = "Genome position",
  ylab = "Copy number ratio",
  main = "Genome CNV profile filtrated"
)

abline(h = 1, lty = 2)
abline(h = 1.5, col = "red", lty = 2)
abline(h = 0.5, col = "blue", lty = 2)




