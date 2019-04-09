# @Date: 8 March 2018
# @Author: Francis Blokzijl
# @Description: Script for generating figures for MutationalPatterns paper

# ----- PACKAGES ----
library("MutationalPatterns")
library("ggplot2")
library("gridExtra")
library("reshape2")
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# ------- DATA ------
# Download the MutPat_object.rds and MutPat_object_invitro.rds from: 
# https://wgs11.op.umcutrecht.nl/mutational_patterns_ASCs/data/

# Specify the directory where the data is downloaded to
data_dir = "..."
# Specify your output directory
out_dir = "..."

# read MutPat objects
MutPat_object = readRDS(paste(data_dir, "MutPat_object.rds", sep = ""))
MutPat_object_invitro = readRDS(paste(data_dir, "MutPat_object_invitro.rds", sep = ""))
# reorder tissue levels for correct plotting order
MutPat_object$tissue = factor(MutPat_object$tissue, levels = c("Colon", "Small Intestine", "Liver"))
MutPat_object_invitro$tissue = factor(MutPat_object_invitro$tissue, levels = c("Colon", "Small Intestine", "Liver"))

# ------ MUTATION MATRICES ------

# mutation matrix
mut_mat_adult = mut_matrix(MutPat_object$vcf, ref_genome)
# mutation matrix with in vitro
mut_mat_invitro = mut_matrix(MutPat_object_invitro$vcf, ref_genome)
# combine mutation matrices of adult and in vitro
mut_mat_plus = cbind(mut_mat_adult, mut_mat_invitro)

# ------- SPECTRUM ------

# get mutation types
type_occurences = mut_type_occurrences(MutPat_object$vcf, ref_genome)
# plot spectrum
plot_spectrum = plot_spectrum(type_occurences, CT = T, by = MutPat_object$tissue)

pdf(paste(out_dir, "spectrum.pdf", sep = ""), width = 7, height = 3, useDingbats = F)
plot_spectrum
dev.off()

# ----- SIGNATURES DE NOVO ------

# Extract 3 signatures from mutation count matrix with in vitro data
nmf_res = extract_signatures(mut_mat_plus, rank = 3)

# Add signature names
colnames(nmf_res$signatures) = c("Signature B", "Signature C", "Signature A")
rownames(nmf_res$contribution) = c("Signature B", "Signature C", "Signature A")

# reorder signatures: A-B-C
order_signatures = c(3,1,2)
nmf_res$signatures = nmf_res$signatures[,order_signatures]
nmf_res$contribution = nmf_res$contribution[order_signatures,]

# plotting
plot_signatures = plot_96_profile(nmf_res$signatures, condensed = T)

pdf(paste(out_dir, "signatures.pdf", sep = ""), width = 7, height = 5, useDingbats = F)
plot_signatures
dev.off()

# Plot signature contribution
# select only adult sample (not in vitro) in order Colon, Liver, Small Intestine
sample_order = c(1:21, 32:45, 22:31)
sample_levels = colnames(mut_mat_plus)[sample_order]
plot_contribution_bar = plot_contribution(nmf_res$contribution[,sample_order], coord_flip = T)

pdf(paste(out_dir, "contribution.pdf", sep = ""), width = 4, height = 8, useDingbats = F)
plot_contribution_bar
dev.off()

# ------ PROOF OF PRINCIPLE NNLS -----

# Signatures are not normalized yet
sigs = nmf_res$signatures
profiles = mut_mat_plus[, sample_order]

fit_res_proof_nnls = fit_to_signatures(profiles, sigs)
# plot contribution barplot
plot_contribution_bar_proof_nnls = plot_contribution(fit_res_proof_nnls$contribution, coord_flip = T)

pdf(paste(out_dir, "proof_nnls_contribution.pdf", sep = ""), width = 10, height = 9, useDingbats = F)
grid.arrange((plot_contribution_bar + ggtitle("NMF")), (plot_contribution_bar_proof_nnls + ggtitle("NNLS")), ncol = 2)
dev.off()

# Calculate correlation
nmf_contribution = nmf_res$contribution[,sample_order]
nnls_contribution = fit_res_proof_nnls$contribution

res_cor = c()
# for all samples
for(i in 1:45)
{
  res_cor[i] = cor(nmf_contribution[,i], nnls_contribution[,i], method = "pearson")
}
mean(res_cor) #0.9787363

# -------- COSMIC CANCER SIGNATURES -------

# download signatures from pan-cancer study Alexandrov et al.
sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
cancer_signatures = read.table(sp_url, sep = "\t", header = T)
# match the order to MP standard in mut_matrix
order = match(row.names(mut_mat_adult), cancer_signatures$Somatic.Mutation.Type)
# reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[order,]
# add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

cos_sim_COSMIC = cos_sim_matrix(cancer_signatures, cancer_signatures)
cosine_heatmap_COSMIC = plot_cosine_heatmap(round(cos_sim_COSMIC,2), plot_values = T, cluster_rows = F)

pdf(paste(out_dir, "cosine_heatmap_COSMIC.pdf", sep = ""), width = 11, height = 10, useDingbats = F)
cosine_heatmap_COSMIC
dev.off()

# similarity of our signatures with cancer cosmic signatures
cos_sim_signatures = cos_sim_matrix(nmf_res$signatures, cancer_signatures)
plot_cos_sim_signatures = plot_cosine_heatmap(cos_sim_signatures, cluster_rows = F, plot_values = T)

pdf(paste(out_dir, "signatures_sim.pdf", sep = ""), width = 12, height = 2, useDingbats = F)
plot_cos_sim_signatures + xlab("COSMIC cancer signatures") + theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()


# ------ FIT TO SIGNATURES -----

# fit mutational profiles to cancer signatures
fit_res = fit_to_signatures(mut_mat_adult, cancer_signatures)

# # timing of algorithm
# # 45 samples
# system.time(fit_to_signatures(mut_mat_adult, cancer_signatures))
# # user  system elapsed 
# # 0.107   0.002   0.110 
# # 90 samples
# system.time(fit_to_signatures(cbind(mut_mat_adult,mut_mat_adult), cancer_signatures))
# # user  system elapsed 
# # 0.179   0.004   0.183
# # 180 samples
# system.time(fit_to_signatures(cbind(mut_mat_adult,mut_mat_adult,mut_mat_adult,mut_mat_adult), cancer_signatures))
# # user  system elapsed 
# # 0.374   0.006   0.382 
# # 360 samples
# mat_180 = cbind(mut_mat_adult,mut_mat_adult,mut_mat_adult,mut_mat_adult)
# system.time(fit_to_signatures(cbind(mat_180, mat_180), cancer_signatures))
# # user  system elapsed 
# # 0.799   0.010   0.814 

# get relative contribution
contribution = t(fit_res$contribution)
# relative contribution
contribution_norm = contribution / rowSums(contribution)
# get maximum contribution over all samples per signature
max_contribution_norm = apply(contribution_norm, 2, function(x) max(x)) 
# get signatures with at least 10% contribution in at least 1 sample
select2 = which(max_contribution_norm > 0.1)
# plot contribution heatmap
plot_fit_contribution_heatmap = plot_contribution_heatmap(fit_res$contribution[select2,sample_order], plot_values = F, cluster_samples = F)

# Fig 2A
pdf(paste(out_dir, "refit_contribution_heatmap.pdf", sep = ""), width = 2.8, height = 6.5, useDingbats = F)
plot_fit_contribution_heatmap + theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()

# cosine similarity matrix between original & reconstructed profiles
cos_sim_ori_rec_cancer = cos_sim_matrix(mut_mat_adult, fit_res$reconstructed)
x = as.matrix(diag(cos_sim_ori_rec_cancer))
colnames(x) = "reconstructed"
plot_fit_reconstructed = plot_cosine_heatmap(x, cluster_rows = F, plot_values = T)
# make as barplot
df = as.data.frame(x)
df$sample = row.names(df)
df$sample = factor(df$sample, levels = df$sample[sample_order])

plot_cosine_ori_rec = ggplot(df, aes(y=reconstructed, x=sample)) + 
  geom_bar(stat="identity", fill = "skyblue4") +
  coord_flip(ylim=c(0.92,1)) +
  ylab("Cosine similarity\n original VS reconstructed") +
  xlab("") +
  # reverse order of the samples such that first is up
  xlim(rev(levels(factor(df$sample)))) +
  theme_bw() +
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank()) +
  geom_hline(aes(yintercept=.95))

# Fig 2B
pdf(paste(out_dir, "cosine_barplot.pdf", sep = ""), width = 2, height = 7, useDingbats = F)
plot_cosine_ori_rec
dev.off()

# get mean 
mean(df$reconstructed) #0.9779207

# take sample 1-a, as this cosine similarity is the lowest
plot_compare_spectra = plot_compare_profiles(mut_mat_adult[,1], fit_res$reconstructed[,1], condensed = T, profile_names = c("Original", "Reconstructed"))

# Fig 2C
pdf(paste(out_dir, "compare_profiles.pdf", sep = ""), width = 7, height = 5.5, useDingbats = F)
plot_compare_spectra
dev.off()

# cluster cancer signatures
hclust_cancer_signatures1 = cluster_signatures(cancer_signatures, method = "average")
cancer_signatures_order1 = colnames(cancer_signatures)[hclust_cancer_signatures1$order]

# cosine similarity with cosmic cancer signatures
cos_sim_samples_signatures = cos_sim_matrix(mut_mat_adult, cancer_signatures)
# plot heatmap
plot_cos_sim_samples_signatures1 = plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cancer_signatures_order1)

# Fig 3
pdf(paste(out_dir,"cosine_heatmap_samples_signatures.pdf", sep = ""), width = 8, height = 6.5, useDingbats = F)
plot_cos_sim_samples_signatures1
dev.off()


# --------- COMPARISON WITH DECONSTRUCTSIGS -------

# install.packages("deconstructSigs")
library(deconstructSigs)

mut_input = as.data.frame(t(mut_mat_adult))
# normalize
mut_input = mut_input / rowSums(mut_input)

# normalize fit result from MutationalPatterns for comparison
fit_res_norm = t(fit_res$contribution)
fit_res_norm = fit_res_norm / rowSums(fit_res_norm)

# for all samples
decon_res = data.frame()  
system.time(
  for(i in 1:ncol(mut_mat_adult))
  {
    res = whichSignatures(tumor.ref = mut_input, 
                          sample.id = colnames(mut_mat_adult)[i], 
                          signatures.ref = signatures.cosmic,
                          signature.cutoff = 0)
    
    decon_res = rbind(decon_res, res$weights)
  })
# user  system elapsed 
# 39.133   5.985  46.339 

# rename signatures to just numeric
colnames(decon_res) = as.character(1:30)

# calculate correlation between the results of the two methods
res_cor = c()
for(i in 1:nrow(decon_res))
{
  res_cor[i] = cor(as.numeric(decon_res[i,]), fit_res_norm[i,])
}
mean(res_cor) #0.98017
median(res_cor) #0.9963108

# check if lapply is faster than for loop
# system.time(
#   lapply(colnames(mut_mat_adult), function(x) whichSignatures(tumor.ref = mut_input, 
#                                                                       sample.id = x, 
#                                                                       signatures.ref = signatures.cosmic,
#                                                                       signature.cutoff = 0)$weights))
# similar timing result
# user  system elapsed 
# 38.045   4.171  42.441

# Plot contribution heatmaps
MP = plot_contribution_heatmap(t(decon_res)[,sample_order], cluster_samples = F) + ggtitle("MutationalPatterns") + theme(axis.text.x = element_text(angle = 0, hjust = 1))
DS = plot_contribution_heatmap(t(fit_res_norm)[,sample_order], cluster_samples = F) + ggtitle("deconstructSigs") + theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Combine in one plot
# Fig S1A
pdf(paste(out_dir, "contribution_MP_DS.pdf", sep = ""), width = 13, height = 7.5, useDingbats = F)
grid.arrange(MP, DS, ncol=2)
dev.off()

# reconstruct mutation count matrix
DS_reconstructed = as.matrix(t(signatures.cosmic)) %*% as.matrix(t(decon_res))
MP_reconstructed = as.matrix(t(signatures.cosmic)) %*% as.matrix(t(fit_res_norm))

# calculate cosine similarity matrix between original & reconstructed
# MutationalPatterns
cos_sim_MP = cos_sim_matrix(mut_mat_adult, MP_reconstructed)
cos_sim_MP = as.matrix(diag(cos_sim_MP))
# DeconstructSigs
cos_sim_DS = cos_sim_matrix(mut_mat_adult, DS_reconstructed)
cos_sim_DS = as.matrix(diag(cos_sim_DS))

# compare cosine similarities between MP & DS
mean(cos_sim_MP) # 0.9779207
mean(cos_sim_DS) # 0.9773019

df = cbind(cos_sim_MP, cos_sim_DS)
colnames(df) = c("MutationalPatterns", "DeconstructSigs")
df = melt(df)
plot_MP_DS_boxplot = ggplot(df, aes(x=Var2, y=value, fill = Var2)) + 
  geom_boxplot() +
  geom_point() +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
  ylim(0.925, 1) +
  ylab("Cosine similarity \noriginal VS reconstructed") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none") 

# Calculate RSS between original and reconstructed for both methods


calculate_RSS = function(profile1, profile2)
{
  # Normalize mutation profiles
  s1_relative = profile1 / sum(profile1)
  s2_relative = profile2 / sum(profile2)
  diff = s1_relative - s2_relative
  
  # residual sum of squares
  RSS = sum(diff^2)
  
  return(RSS)
}

RSS_MP = c()
for(i in 1:45)
{
  new_RSS = calculate_RSS(mut_mat_adult[,i], MP_reconstructed[,i])
  RSS_MP = c(RSS_MP, new_RSS)
}

RSS_DS = c()
for(i in 1:45)
{
  new_RSS = calculate_RSS(mut_mat_adult[,i], DS_reconstructed[,i])
  RSS_DS = c(RSS_DS, new_RSS)
}

# compare RSS between MP & DS
mean_RSS_MP = mean(RSS_MP) # 0.001379043 
mean_RSS_DS = mean(RSS_DS) # 0.001400141

# Format
format(mean_RSS_MP, scientific = TRUE, digits = 3) # "1.38e-03"
format(mean_RSS_DS, scientific = TRUE, digits = 3) # "1.4e-03"

# Fig S1B
pdf(paste(out_dir,"boxplot_MP_DS.pdf", sep = ""), width = 3, height = 4, useDingbats = F)
plot_MP_DS_boxplot
dev.off()

# ------- TRANSCRIPTIONAL STRAND BIAS ------

# get known genes table from UCSC for hg19 using
# source("https://bioconductor.org/biocLite.R")
# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
# reduce
genes_hg19 = reduce(genes_hg19)
# make stranded mutation matrix
mut_mat_s_adult = mut_matrix_stranded(MutPat_object$vcf, ref_genome, genes_hg19)
# in vitro
mut_mat_s_invitro = mut_matrix_stranded(MutPat_object_invitro$vcf, ref_genome, genes_hg19)
# combiend
mut_mat_s_plus = cbind(mut_mat_s_adult, mut_mat_s_invitro)
# count
strand_counts = strand_occurrences(mut_mat_s_adult, by = MutPat_object$tissue)
# plot
plot_strand_assym = plot_strand(strand_counts, mode = "relative")

strand_bias = strand_bias_test(strand_counts)
plot_strand_bias = plot_strand_bias(strand_bias)

# Fig 4
pdf(paste(out_dir, "strand_bias.pdf", sep = ""), width = 8, height = 6, useDingbats = F)
grid.arrange(plot_strand_assym, plot_strand_bias)
dev.off()

# ------ SIGNATURES WITH STRAND BIAS -----

nmf_res_strand = extract_signatures(mut_mat_s_plus, rank = 3)
plot_192_profile(nmf_res_strand$signatures)
colnames(nmf_res_strand$signatures) = c("Signature A", "Signature B", "Signature C")

plot_signatures_strand2 = plot_192_profile(nmf_res_strand$signatures, condensed = T)
plot_signatures_strand_bias = plot_signature_strand_bias(nmf_res_strand$signatures)

# Fig 4
pdf(paste(out_dir,"signatures_strand_bias.pdf", sep = ""), width = 11, height = 5, useDingbats = F)
grid.arrange(plot_signatures_strand2, plot_signatures_strand_bias, ncol=2, widths=c(5,2))
dev.off()

# ------- REPLICATION STRAND BIAS -------

repli_file = system.file("extdata/ReplicationDirectionRegions.bed",package = "MutationalPatterns")
repli_strand = read.table(repli_file, header = TRUE)
# Store in GRanges object
repli_strand_granges = GRanges(seqnames = repli_strand$Chr,
                               ranges = IRanges(start = repli_strand$Start + 1,end = repli_strand$Stop),
                               strand_info = repli_strand$Class)


# UCSC seqlevelsstyle
seqlevelsStyle(repli_strand_granges) = "UCSC"
repli_strand_granges

# Make mutation count matrix with transcriptional strand information 
# (96 trinucleotides * 2 strands = 192 features).
mut_mat_s_rep <- mut_matrix_stranded(MutPat_object$vcf, ref_genome, repli_strand_granges,mode = "replication")
# View
mut_mat_s_rep[1:5, 1:5]
# Count the number of mutations on each strand, per tissue, per mutation type
strand_counts_rep <- strand_occurrences(mut_mat_s_rep, by=MutPat_object$tissue)
head(strand_counts_rep)
# Perform Poisson test for strand asymmetry significance testing:
strand_bias_rep <- strand_bias_test(strand_counts_rep)
strand_bias_rep
# Plot the mutation spectrum with strand distinction:
ps1 <- plot_strand(strand_counts_rep, mode = "relative")
#Plot the effect size (log2(untranscribed/transcribed) of the strand bias. Asteriks indicate significant strand bias.
ps2 <- plot_strand_bias(strand_bias_rep)
# Combine the plots into one figure:
grid.arrange(ps1, ps2)

pdf(paste(out_dir, "rep_strand_bias.pdf", sep = ""), width = 8, height = 7, useDingbats = F)
grid.arrange(ps1, ps2)
dev.off()

# ------ RAINFALL PLOT ------

chromosomes = names(genome(vcf_adult[[1]])[1:22])

pdf(paste(out_dir, "rainfalls.pdf", sep = ""), useDingbats = F)
for(i in 1:45){
  plot(plot_rainfall(MutPat_object$vcf[[i]], title = MutPat_object$sample_id[i], chromosomes = chromosomes))}
dev.off()

# choose interesting one: sample 14b
pdf(paste(out_dir, "rainfall_14b.pdf", sep = ""), width = 11, height = 4, useDingbats = F)
plot(plot_rainfall(MutPat_object$vcf[[23]], chromosomes = chromosomes, cex=2))
dev.off()

# ------ GENOMIC DISTRIBUTION ------

# download promoter annotations with biomaRt
# this data is included in the package
promoter_g = readRDS(system.file("states/promoter_g_data.rds", package="MutationalPatterns"))
seqlevelsStyle(promoter_g) =  "UCSC"
# download genes and non-genic genomic regions
hg19_Granges = GRanges(seqnames = seqnames(get(ref_genome))[1:22], IRanges(start=rep(1,22), end=seqlengths(get(ref_genome))[1:22]))
hg19_nongenic = setdiff(hg19_Granges, genes_hg19, ignore.strand=T)
# Combine all genomic regions (GRanges objects) in a named list:
regions = list(promoter_g, genes_hg19, hg19_nongenic)
names(regions) = c("Promoter", "Genes", "Non-genic")

# genomic distribution testing
distr = genomic_distribution(MutPat_object$vcf, MutPat_object$surveyed, regions)
distr_test = enrichment_depletion_test(distr, by = MutPat_object$tissue)
# change order for plotting
distr_test$by = factor(distr_test$by, levels = c("Colon", "Small Intestine", "Liver"))

# Fig 4
pdf(paste(out_dir,"genomic_distribution.pdf", sep = ""), width = 7, height = 6, useDingbats = F)
plot_enrichment_depletion(distr_test)
dev.off()
