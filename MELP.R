################################################################################

####################### Gene Expression Analysis ###############################

################################################################################
library(ggplot2)
library(tidyverse)
library(coRdon)
library(ggpubr)


# Import sequence data & generate codon tables

################################################################################

####################################### LPAI ###################################

################################################################################

# Low Pathogenic Avian Influenza HA
LPAI_HA_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\chicken_LPAI_HA_733.fna"
)
LPAI_HA_ct <- codonTable(LPAI_HA_dna)

# Low Pathogenic Avian Influenza NP
LPAI_NP_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\chicken_LPAI_NP.fna"
)
LPAI_NP_ct <- codonTable(LPAI_NP_dna)

# Low Pathogenic Avian Influenza M1 (Host = Chicken)
LPAI_M1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\chicken_LPAI_MP.fasta"
)
LPAI_M1_ct <- codonTable(LPAI_M1_dna)

# Low Pathogenic Avian Influenza NS1 (Host = Duck)
LPAI_NS1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\duck_LPAI_NS1.txt"
)
LPAI_NS1_ct <- codonTable(LPAI_NS1_dna)

################################################################################

#################################### HPAI ######################################

################################################################################

# High Pathogenic Avian Influenza HA (Host = Chicken)
HPAI_HA_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\chicken_HPAI_HA_8684.fna"
)
HPAI_HA_ct <- codonTable(HPAI_HA_dna)

# High Pathogenic Avian Influenza NP (Host = Chicken)
HPAI_NP_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\chicken_HPAI_NP.fna"
)
HPAI_NP_ct <- codonTable(HPAI_NP_dna)

# High Pathogenic Avian Influenza M1 (Host = Chicken)
HPAI_M1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\chicken_HPAI_MP.fasta"
)
HPAI_M1_ct <- codonTable(HPAI_M1_dna)

# High Pathogenic Avian Influenza NS1 (Host = Duck)
HPAI_NS1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\duck_HPAI_NS1.txt"
)
HPAI_NS1_ct <- codonTable(HPAI_NS1_dna)

# Reference Sequences
H5N1_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\HA\\hemagglutinin.fasta"
)
H5N1_ha <- codonTable(H5N1_ha_dna)

H5N1_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\NP\\nucleocapsid.fasta"
)
H5N1_np <- codonTable(H5N1_np_dna)

H5N1_m1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\MP1\\matrix protein 1.fasta"
)
H5N1_m1 <- codonTable(H5N1_m1_dna)

H5N1_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\NSP1\\nonstructural protein 1.fasta"
)
H5N1_ns1 <- codonTable(H5N1_ns1_dna)


# Run MELP analysis to determine difference in HA, NP, or M1 gene expression
melp_ha <- MELP(H5N1_ha, subsets = list(LPAI_HA_ct, HPAI_HA_ct))
melp_np <- MELP(H5N1_np, subsets = list(LPAI_NP_ct, HPAI_NP_ct))
melp_m1 <- MELP(H5N1_m1, subsets = list(LPAI_M1_ct, HPAI_M1_ct))
melp_ns1 <- MELP(H5N1_ns1, subsets = list(LPAI_NS1_ct, HPAI_NS1_ct))

# Convert the data to a data frame & change column names
melp_ha_df <- as.data.frame(melp_ha)
colnames(melp_ha_df) <- c("LPAI", "HPAI")
melp_ha_df <- melp_ha_df %>% mutate(Gene = "HA")
melp_ha_df <- melp_ha_df %>% 
  pivot_longer(c("LPAI", "HPAI"), names_to = "Pathogenicity", values_to = "Predicted Expression")

melp_np_df <- as.data.frame(melp_np)
colnames(melp_np_df) <- c("LPAI", "HPAI")
melp_np_df <- melp_np_df %>% mutate(Gene = "NP")
melp_np_df <- melp_np_df %>% 
  pivot_longer(c("LPAI", "HPAI"), names_to = "Pathogenicity", values_to = "Predicted Expression")

melp_m1_df <- as.data.frame(melp_m1)
colnames(melp_m1_df) <- c("LPAI", "HPAI")
melp_m1_df <- melp_m1_df %>% mutate(Gene = "M1")
melp_m1_df <- melp_m1_df %>% 
  pivot_longer(c("LPAI", "HPAI"), names_to = "Pathogenicity", values_to = "Predicted Expression")

melp_ns1_df <- as.data.frame(melp_ns1)
colnames(melp_ns1_df) <- c("LPAI", "HPAI")
melp_ns1_df <- melp_ns1_df %>% mutate(Gene = "NS1")
melp_ns1_df <- melp_ns1_df %>% 
  pivot_longer(c("LPAI", "HPAI"), names_to = "Pathogenicity", values_to = "Predicted Expression")
melp_ns1_df <- melp_ns1_df[melp_ns1_df$`Predicted Expression` > 0,]
# Not working properly, LPAI data not showing expression data. Excluding for now.


melp_df <- rbind(melp_ha_df, melp_np_df, melp_m1_df)

# Plot and visualize the data
ggplot(data = melp_df, aes(x = Gene, y = `Predicted Expression`, fill = Pathogenicity)) +
  geom_boxplot() +
  labs(title = "MELP Analysis of HPAI vs. LPAI Strains") +
  stat_compare_means(method = "t.test")

# Now add in other hosts to see if expression varies by host
# HA
ha_human_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\human_HPAI_HA_60.fna"
)
ha_human_ct <- codonTable(ha_human_dna)

ha_cow_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\dairycow_HPAI_HA_3008.fna"
)
ha_cow_ct <- codonTable(ha_cow_dna)

#NP
np_human_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\human_np.fna"
)
np_human_ct <- codonTable(np_human_dna)

np_cow_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\cow_np.fna"
)
np_cow_ct <- codonTable(np_cow_dna)

# M1
m1_human_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\Human_HPAI_MP.fna"
)
m1_human_ct <- codonTable(m1_human_dna)

m1_cow_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\dairycow_HPAI_MP.fna"
)
m1_cow_ct <- codonTable(m1_cow_dna)

# NS1
ns1_human_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\human_HPAI_ns1.fasta"
)
ns1_human_ct <- codonTable(ns1_human_dna)

ns1_cow_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\cow_ns1.fasta"
)
ns1_cow_ct <- codonTable(ns1_cow_dna)
ns1_chicken_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\chicken_ns1.fasta"
)
ns1_chicken_ct <- codonTable(ns1_chicken_dna)

# HA MELP by Host
ha_by_host <- MELP(H5N1_ha, subsets = list(HPAI_HA_ct, ha_cow_ct, ha_human_ct))
ha_by_host <- as.data.frame(ha_by_host)
colnames(ha_by_host) <- c("Chicken", "Cow", "Human")
ha_by_host <- ha_by_host %>%
  pivot_longer(cols = c("Chicken", "Cow", "Human"), names_to = "Host", values_to = "Predicted Expression")
ha_by_host <- ha_by_host %>% mutate(Gene = "HA")
ha_melp_mean <- mean(ha_by_host$`Predicted Expression`)
# HA mean = 0.890

# NP MELP by Host
np_by_host <- MELP(H5N1_np, subsets = list(HPAI_NP_ct, np_cow_ct, np_human_ct))
np_by_host <- as.data.frame(np_by_host)
colnames(np_by_host) <- c("Chicken", "Cow", "Human")
np_by_host <- np_by_host %>%
  pivot_longer(cols = c("Chicken", "Cow", "Human"), names_to = "Host", values_to = "Predicted Expression")
np_by_host <- np_by_host %>% mutate(Gene = "NP")
np_melp_mean <- mean(np_by_host$`Predicted Expression`)
# NP mean = 0.842

# M1 MELP by Host
m1_by_host <- MELP(H5N1_m1, subsets = list(HPAI_M1_ct, m1_cow_ct, m1_human_ct))
m1_by_host <- as.data.frame(m1_by_host)
colnames(m1_by_host) <- c("Chicken", "Cow", "Human")
m1_by_host <- m1_by_host %>%
  pivot_longer(cols = c("Chicken", "Cow", "Human"), names_to = "Host", values_to = "Predicted Expression")
m1_by_host <- m1_by_host %>% mutate(Gene = "M1")
m1_melp_mean <- mean(m1_by_host$`Predicted Expression`)
# M1 mean = 0.657

# NS1 MELP by Host
ns1_by_host <- MELP(H5N1_ns1, subsets = list(HPAI_NS1_ct, ns1_cow_ct, ns1_human_ct))
ns1_by_host <- as.data.frame(ns1_by_host)
colnames(ns1_by_host) <- c("Chicken", "Cow", "Human")
ns1_by_host <- ns1_by_host %>%
  pivot_longer(cols = c("Chicken", "Cow", "Human"), names_to = "Host", values_to = "Predicted Expression")
ns1_by_host <- ns1_by_host %>% mutate(Gene = "NS1")
ns1_melp_mean <- mean(ns1_by_host$`Predicted Expression`)
# Not working properly... excluding this portion of the analysis

melp_by_host <- rbind(ha_by_host, np_by_host, m1_by_host, ns1_by_host)
melp_by_host <- melp_by_host[melp_by_host$`Predicted Expression` > 0.3,]

ggplot(data = melp_by_host, aes(x = Host, y = `Predicted Expression`, fill = Gene)) +
  geom_boxplot() +
  labs(title = "Predicted Expression of HA, NP, NS1, and M1 by Host")

################################################################################

################### MELP of Highly Divergent Sequences #########################

################################################################################

# Import sequences which fell further away on the tree
# NOTE - Sequence selection was ambiguous with no particular pattern or parameters.
# Sequences are derived from human hosts

mutated_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Analysis\\Alignments\\Mutated_NP_seqs.fasta"
)
mutated_np_ct <- codonTable(mutated_np_dna)

mutated_m1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Analysis\\Alignments\\Mutated_M1_seqs.fasta"
)
mutated_m1_ct <- codonTable(mutated_m1_dna)

mutated_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Analysis\\Alignments\\mutated_ha.fasta"
)
mutated_ha_ct <- codonTable(mutated_ha_dna)

mutated_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Analysis\\Alignments\\Mutated_NS1_seqs.fasta"
)
mutated_ns1_ct <- codonTable(mutated_ns1_dna)

# Run MELP with comparison to human host sequences
mutated_np_melp <- MELP(np_human_ct, subsets = list(mutated_np_ct))
mutated_m1_melp <- MELP(m1_human_ct, subsets = list(mutated_m1_ct))
mutated_ha_melp <- MELP(ha_human_ct, subsets = list(mutated_ha_ct))
mutated_ns1_melp <- MELP(ns1_human_ct, subsets = list(mutated_ns1_ct))

# Restructure the data
mutated_np_df <- as.data.frame(mutated_np_melp)
mutated_np_df <- mutated_np_df %>% mutate(Gene = "NP")
colnames(mutated_np_df) <- c("Predicted Expression", "Gene")

mutated_m1_df <- as.data.frame(mutated_m1_melp)
mutated_m1_df <- mutated_m1_df %>% mutate(Gene = "M1")
colnames(mutated_m1_df) <- c("Predicted Expression", "Gene")

mutated_ha_df <- as.data.frame(mutated_ha_melp)
mutated_ha_df <- mutated_ha_df %>% mutate(Gene = "HA")
colnames(mutated_ha_df) <- c("Predicted Expression", "Gene")

mutated_ns1_df <- as.data.frame(mutated_ns1_melp)
mutated_ns1_df <- mutated_ns1_df %>% mutate(Gene = "NS1")
colnames(mutated_ns1_df) <- c("Predicted Expression", "Gene")

mutated_genes <- rbind(mutated_np_df, mutated_m1_df, mutated_ha_df, mutated_ns1_df)
mutated_genes <- mutated_genes %>% mutate(Origin_Set = "Mutated")
mutated_genes <- mutated_genes[mutated_genes$`Predicted Expression` > 0,]

# Make a new dataframe for the plot
human_melp <- melp_by_host[melp_by_host$Host == "Human",]
human_melp <- subset(human_melp, select = -Host)
human_melp <- human_melp %>% mutate(Origin_Set = "All")

overall_v_mutated <- rbind(human_melp, mutated_genes)

# Plot the original MELP values against the ones for the samples which were more distant on the phylogenic tree
ggplot(data = overall_v_mutated, aes(x = Gene, y = `Predicted Expression`, fill = Origin_Set)) +
  geom_boxplot() +
  labs(title = "Comparison of Total MELP and Mutated MELP Values") +
  stat_compare_means(method = "t.test", label.x = 0.5, label.y = 1.15)

