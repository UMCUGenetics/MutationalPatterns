# USER MANUAL
library(MutationalPatterns)
library(gridExtra)
vcf_files = list.files(system.file("extdata", package="MutationalPatterns"), full.names = T)
sample_names = c("colon1", "colon2", "colon3", "intestine1", "intestine2", "intestine3", "liver1", "liver2", "liver3")
vcfs = read_vcf(vcf_files, sample_names)

ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
# metadata
tissue = c("colon", "colon", "colon", "intestine", "intestine", "intestine", "liver", "liver", "liver")

# ------ BASIC FUNCTIONS ----

muts = get_muts(vcfs[[1]])
types = get_types(vcfs[[1]])
mut_context = get_mut_context(vcfs[[1]], ref_genome)
type_context = get_type_context(vcfs[[1]], ref_genome)
type_occurences = mut_type_occurences(vcfs[1], ref_genome)
type_occurences = mut_type_occurences(vcfs, ref_genome)


# ----- SPECTRA ------

# plot spectrum over all samples
plot1 = plot_spectrum(type_occurences)
# with distinction of C>T at CpG sites
plot2 = plot_spectrum(type_occurences, CT = T)
# without legend
plot3 = plot_spectrum(type_occurences, CT = T, legend = F)
# define own colors for plotting
my_colors = c("pink", "orange", "blue", "lightblue", "green", "red", "purple")
plot4 = plot_spectrum(type_occurences, CT = T, legend = T, colors = my_colors)
# plot spectrum per tissue
plot5 = plot_spectrum(type_occurences, by = tissue, CT = T)

spectra1 = grid.arrange(plot1, plot2, plot3, ncol=3, widths=c(3,3,1.8))
spectra2 = grid.arrange(plot5, plot4, ncol=2, widths=c(3,2))



# ------ 96 SPECTRUM -----

# make 96 count matrix
mut_matrix = make_mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
# plot 96 profile of three samples
plot = plot_96_profile(mut_matrix[,c(1,4,7)]) 
plot

# ------ COMPARE PROFLIES ------

# get a p-value for the significance of a difference between two profiles
comparison <- profile_bootstrap_comparison(mut_matrix[,1], mut_matrix[,4], random.seed = 123)
comparison$overallPvalue
plot_96_profile(mut_matrix[,c(1,4)]) 

# get the mutation types that significantly differ between two profiles
comparison <- profile_bootstrap_comparison(mut_matrix[,1], mut_matrix[,7], random.seed = 123)
comparison$mutTypePvalues_corrected[comparison$mutTypePvalues_corrected < 0.05]
plot_96_profile(mut_matrix[,c(1,7)]) 


# ------ EXTRACT SIGNATURES ------

# estimate rank
estimate_rank(mut_matrix, rank_range = 2:5, nrun = 50)
# extract 3 signatures
nmf_res = extract_signatures(mut_matrix, rank = 3)
# provide signature names (optional)
colnames(nmf_res$signatures) = c("Signature A", "Signature B" , "Signature C")
# plot signatures
plot_96_profile(nmf_res$signatures)

# provide signature names (optional)
rownames(nmf_res$contribution) = c("Signature A", "Signature B" , "Signature C")
# plot signature contribution
plot_contribution(nmf_res$contribution, coord_flip = T)

# plot reconstructed mutation profile
plot_96_profile(nmf_res$reconstructed)

# compare reconstructed 96 profile of sample 1 with orignal profile
plot_compare_profiles(mut_matrix[,1], nmf_res$reconstructed[,1], profile_names = c("Original", "Reconstructed"))

# control if the reconstruction significantly matches the original
profile_bootstrap_comparison(mut_matrix[,1], nmf_res$reconstructed[,1])$overallPvalue


# ------ REFIT SIGNATURES ------

# download signatures from pan-cancer study Alexandrov et al.
url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
cancer_signatures = read.table(url, sep = "\t", header = T)
# reorder (to make the order of the trinucleotide changes the same)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
# only signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
# plot signature 1
plot_96_profile(cancer_signatures[,1,drop=F], ymax = 0.2)

fit_res = fit_to_signatures(mut_matrix, cancer_signatures)
# select signatures with some contribution
select = which(rowSums(fit_res$contribution) > 0)
# plot contribution
plot_contribution(fit_res$contribution[select,], coord_flip = T)

# compare reconstructed from refit with original profile
plot_compare_profiles(mut_matrix[,1], fit_res$reconstructed[,1], profile_names = c("Original", "Reconstructed \n cancer signatures"))

# control if the reconstruction significantly matches the original
profile_bootstrap_comparison(mut_matrix[,1], fit_res$reconstructed[,1])$overallPvalue



# ------ RAINFALL PLOT ------

# define chromosomes of interest
chromosomes = seqnames(get(ref_genome))[1:22]
# along chromosome
plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes, cex = 1.5)
# for chromosome 1
plot_rainfall(vcfs[[1]], title = names(vcfs[1]), ref_genome = ref_genome, chromosomes = chromosomes[1], cex = 2)


# ------ GENOMIC DISTRIBUTION -------

bed_files = list.files("~/Documents/Organoids/ASC_multi_tissues/data/genomic_regions/", full.names = T)
region_names = c("H3k27ac", "H3K9me3", "early", "intermediate", "late", "exonic")

regions = bed_to_granges(bed_files, region_names = region_names)

surveyed_files = list.files("~/Documents/Organoids/ASC_multi_tissues/data/surveyed/colon/", full.names = T)
surveyed_consensus = bed_to_granges(surveyed_files[1], "surveyed_consensus")
surveyed_list = c(surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus, surveyed_consensus)

# GENOMIC DISTRIBUTION
x = genomic_distribution(vcfs, surveyed_list, regions)
df = enrichment_depletion_test(x)
# plotting is incorrect when significance asterisks is missing
select = which(!(df$significant == "*"))
plot_enrichment_depletion(df)