library(rtracklayer)
library(Gviz)
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
library(ggplot2)

# Importa o vcf
vcf <- readVcf("/home/janicerpaula/Bioinformática/TreinamentoAnaliseDNA/HCMV/vcf_results/variants.vcf")
# Converte vcf para GRanges
vcf_gr <- rowRanges(vcf)
# Confere nome do cromossomo no vcf
seqlevels(vcf_gr)

# Importa o arquivo gff
gff <- import("~/Bioinformática/TreinamentoAnaliseDNA/HCMV/SpnEff/sequence.gff3")
# Verifica se foi gerado corretamente, Aqui, deve ser gerado um GRanges - Genomic Ranges
class(gff)
# Verifica os levels
unique(gff$type)
# Filtra os genes, CDS, mRNAs
genes <- gff[gff$type == "gene"]
cds <- gff[gff$type == "CDS"]
mrna <- gff[gff$type == "mRNA"]
# Verifica o seqname
seqlevels(gff)

# Nomes diferentes, precisam estar com memso nome do cromossomo
seqlevels(vcf_gr) <- seqlevels(gff)

# Verificar posição das variantes no gff
range(gff)
# Verificar posição das variantes no vcf
range(vcf_gr)

# Ver a primeira variante
vcf_gr[1]

# Ver o primeiro gene
genes[1]

# Teste Manual
subsetByOverlaps(vcf_gr, genes)


# Encontra as variantes dentro de genes
hits <- findOverlaps(vcf_gr, genes)


# Extrair resultados
variants_in_genes <- data.frame(
  variant = queryHits(hits),
  gene = genes$Name[subjectHits(hits)]
)

# Ver o resultado
head(variants_in_genes)

# Contar variantes por gene
table(variants_in_genes$gene)

# Tabela final para o relatório
variants_summary <- variants_in_genes %>%
  count(gene, sort = TRUE)

# Gerar plot top genes mutados
top_genes <- head(variants_summary,20)

ggplot(top_genes, aes(x=reorder(gene,n), y=n)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  labs(
    x="Gene",
    y="Number of Variants",
    title="Genes with highest number of variants"
  )

# Distribuição das variantes no genoma (com Gviz)
options(ucscChromosomeNames = FALSE)

# GeneTrack
geneTrack <- GeneRegionTrack(
  genes,
  genome = "HCMV",
  chromosome = "NC_002978.6",
  name = "Genes",
  transcriptAnnotation = "symbol",
  stacking = "dense",
  col = "darkblue",
  fill = "lightblue"
)

# Gerar o Variant Track Annotation
variantTrack <- AnnotationTrack(
  vcf_gr,
  genome = "HCMV",
  chromosome = "NC_002978.6",
  name = "Variants",
  stacking = "dense",
  shape = "box",
  col = "red",
  fill = "red"
)

# Densidade das mutações
variantDensity <- DataTrack(
  start = start(vcf_gr),
  end = end(vcf_gr),
  data = rep(1, length(vcf_gr)),
  genome = "HCMV",
  chromosome = "NC_002978.6",
  type = "histogram",
  name = "Variant Density"
)

# Track do eixo genomico
axisTrack <- GenomeAxisTrack()

# gera o gviz
geneTrack <- GeneRegionTrack(
  genes,
  genome = "HCMV",
  chromosome = "NC_002978.6",
  name = "Genes"
)

# Plot final - muito grande não tem como rodar
plotTracks(
  list(axisTrack, geneTrack, variantTrack),
  from = 1,
  to = max(end(genes)),
  transcriptAnnotation = "symbol"
)

# Plot das densidades das mutações, apenas
plotTracks(
  list(axisTrack, geneTrack, variantDensity),
  from = 1,
  to = max(end(genes))
)


# NORMALIZAR POR TAMANHO DO GENE
# Calcular tamanho dos genes
gene_lengths <- width(genes)

gene_length_df <- data.frame(
  gene = genes$gene,
  length = gene_lengths
)

# Juntar com variants_summary
variants_density <- variants_summary %>%
  left_join(gene_length_df, by = "gene") %>%
  mutate(mutation_density = n / length)

# Top genes por densidade
variants_density <- variants_density %>%
  arrange(desc(mutation_density))

# Plot
top_density <- head(variants_density, 20)

ggplot(top_density, aes(x = reorder(gene, mutation_density), y = mutation_density)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Gene",
    y = "Mutation Density",
    title = "Genes with highest mutation density"
  )




