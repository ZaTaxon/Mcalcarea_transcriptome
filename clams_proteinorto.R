library(stringr)
library(FactoMineR)

# read file output from proteinortho with all 5 species.
ort <- read.delim("clams.proteinortho.tsv", header = TRUE, na.strings = "*")

# extract orthologies that are in all species
ort_all <- ort[complete.cases(ort),]

# count number of ortholog copies per species
Ngene_arctica <- str_count(ort_all$Aislandica_prot.fa, pattern = ',')+1
N_ortogenes <- data.frame(Aislandica = Ngene_arctica)
N_ortogenes$Harctica <- str_count(ort_all$Harctica_prot.fa, pattern = ',')+1
N_ortogenes$Marenaria <- str_count(ort_all$Marenaria_prot.fa, pattern = ',')+1
N_ortogenes$Lbalthica <- str_count(ort_all$Mbalthica_prot.fa, pattern = ',')+1
N_ortogenes$Mcalcarea <- str_count(ort_all$Mcalcarea_prot.fa, pattern = ',')+1

# read file with total protein sequencies (being input from proteinortho)
N_prot <- read.csv("prot_all_seq.csv", header = T)

# normalize orthology numbers to total protein number per species
N_ortogenes_weight <- N_ortogenes / t(N_prot$N_prot)

# transpose df with orthology numbers for PCA
N_ortogenes_w_t <- t(N_ortogenes_weight)

# PCA by normalized numbers of orthologies
res.mpca_w <- PCA(N_ortogenes_w_t, scale.unit = T, ncp=10, graph=T)

pdf(file = "PCA_orthologies.pdf")
plot(res.mpca_w, cex = 2, cex.axis = 2, sex.lab = 3)
dev.off()

# correlation orthologies with principal components
res.dimdesc <- dimdesc(res.mpca_w, axes = c(1,2))

# save component loadings to files
write.csv(res.dimdesc$Dim.1$quanti, file = "orto_dim1_load.csv")
write.csv(res.dimdesc$Dim.2$quanti, file = "orto_dim2_load.csv")

# count mean orthology number in subtidal and tidal groups (depth groups)
str(N_ortogenes_weight)
N_ortogenes_weight$subtidal_mean <- rowMeans(N_ortogenes_weight[,c(1,2,5)])
N_ortogenes_weight$tidal_mean <- rowMeans(N_ortogenes_weight[,c(3,4)])

#count ratio of orthologies numbers in depth groupth 
N_ortogenes_weight$ratio_ts <- N_ortogenes_weight$tidal_mean / N_ortogenes_weight$subtidal_mean
N_ortogenes_weight$ratio_st <- N_ortogenes_weight$subtidal_mean / N_ortogenes_weight$tidal_mean

# plot ratio tidal to subtidal
pdf(file = "tidal_subtidal_ratio.pdf")
plot((N_ortogenes_weight$ratio_ts), xlab = "ortogroup number", ylab = "tidal to subtidal ratio")
abline(h = 5, col = "red")
dev.off()

range(N_ortogenes_weight$ratio_ts)

# plot ratio subtidal to tidal 
pdf(file = "subtidal_tidal_ratio.pdf")
plot((N_ortogenes_weight$ratio_st), xlab = "ortogroup number", ylab = "subtidal to tidal ratio")
abline(h = 5, col = "red")
dev.off()

range(N_ortogenes_weight$ratio_st)

# extract orthologies with high ratio (> 5) subtidal to tidal
subt_max <- ort_all[N_ortogenes_weight$ratio_st >5,]
subt_max$ratio_st <- N_ortogenes_weight$ratio_st[N_ortogenes_weight$ratio_st >5]

# save orthologies with high ratio (> 5) subtidal to tidal
write.csv(subt_max, file = "subtidal_max.csv")


# see numbers of orthology copies per species for orthologies with high ratio
N_ortogenes_weight[N_ortogenes_weight$ratio_st >5,]*100000

# extract orthologies with high ratio (> 5)  tidal to subtidal
tid_max <- ort_all[N_ortogenes_weight$ratio_ts >5,]
tid_max$ratio_ts <- N_ortogenes_weight$ratio_ts[N_ortogenes_weight$ratio_ts >5]
# save orthologies with high ratio (> 5)  tidal to subtidal
write.csv(tid_max, file = "tidal_max.csv")