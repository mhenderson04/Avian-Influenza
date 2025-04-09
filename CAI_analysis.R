library(ggplot2)
library(tidyverse)
library(coRdon)

# Import chicken and human sequences to use as reference for CAI
human_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\human_genes_from_ensembl.fasta"
)
human_ct <- codonTable(human_dna)

chicken_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\Gallus_gallus.fasta"
)
chicken_ct <- codonTable(chicken_dna)

################################################################################

#################################### H5N1 ######################################

################################################################################

# Import sequence files:
H5N1_pb2_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\PB2\\PB2.fasta"
)
H5N1_pb1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\PB1\\PB1.fasta"
)
H5N1_pa_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\PA\\PA.fasta"
)
H5N1_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\HA\\hemagglutinin.fasta"
)
H5N1_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\NP\\nucleocapsid.fasta"
)
H5N1_na_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\NA\\neuraminidase.fasta"
)
H5N1_m1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\MP1\\matrix protein 1.fasta"
)
H5N1_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N1_genes\\NSP1\\nonstructural protein 1.fasta"
)

# Convert to codon tables
H5N1_pb2 <- codonTable(H5N1_pb2_dna)
H5N1_pb1 <- codonTable(H5N1_pb1_dna)
H5N1_pa <- codonTable(H5N1_pa_dna)
H5N1_ha <- codonTable(H5N1_ha_dna)
H5N1_np <- codonTable(H5N1_np_dna)
H5N1_na <- codonTable(H5N1_na_dna)
H5N1_m1 <- codonTable(H5N1_m1_dna)
H5N1_ns1 <- codonTable(H5N1_ns1_dna)

# Establish column names for the data frame
rnames <- c("CAI-Chicken", "CAI-Human")

# Run the CAI analyses, convert to dataframes, manipulate df columns
h5n1_pb2_cai <- CAI(H5N1_pb2, subsets = list(chicken_ct, human_ct))
h5n1_pb2_df <- as.data.frame(h5n1_pb2_cai)
colnames(h5n1_pb2_df) <- rnames
h5n1_pb2_df <- h5n1_pb2_df %>% mutate(Gene = "PB2", Subtype = "H5N1") #Add columns for the boxplot later

h5n1_pb1_cai <- CAI(H5N1_pb1, subsets = list(chicken_ct, human_ct))
h5n1_pb1_df <- as.data.frame(h5n1_pb1_cai)
colnames(h5n1_pb1_df) <- rnames
h5n1_pb1_df <- h5n1_pb1_df %>% mutate(Gene = "PB1", Subtype = "H5N1")

h5n1_pa_cai <- CAI(H5N1_pa, subsets = list(chicken_ct, human_ct))
h5n1_pa_df <- as.data.frame(h5n1_pa_cai)
colnames(h5n1_pa_df) <- rnames
h5n1_pa_df <- h5n1_pa_df %>% mutate(Gene = "PA", Subtype = "H5N1")

h5n1_ha_cai <- CAI(H5N1_ha, subsets = list(chicken_ct, human_ct))
h5n1_ha_df <- as.data.frame(h5n1_ha_cai)
colnames(h5n1_ha_df) <- rnames
h5n1_ha_df <- h5n1_ha_df %>% mutate(Gene = "HA", Subtype = "H5N1")

h5n1_np_cai <- CAI(H5N1_np, subsets = list(chicken_ct, human_ct))
h5n1_np_df <- as.data.frame(h5n1_np_cai)
colnames(h5n1_np_df) <- rnames
h5n1_np_df <- h5n1_np_df %>% mutate(Gene = "NP", Subtype = "H5N1")

h5n1_na_cai <- CAI(H5N1_na, subsets = list(chicken_ct, human_ct))
h5n1_na_df <- as.data.frame(h5n1_na_cai)
colnames(h5n1_na_df) <- rnames
h5n1_na_df <- h5n1_na_df %>% mutate(Gene = "NA", Subtype = "H5N1")

h5n1_m1_cai <- CAI(H5N1_m1, subsets = list(chicken_ct, human_ct))
h5n1_m1_df <- as.data.frame(h5n1_m1_cai)
colnames(h5n1_m1_df) <- rnames
h5n1_m1_df <- h5n1_m1_df %>% mutate(Gene = "M1", Subtype = "H5N1")

h5n1_ns1_cai <- CAI(H5N1_ns1, subsets = list(chicken_ct, human_ct))
h5n1_ns1_df <- as.data.frame(h5n1_ns1_cai)
colnames(h5n1_ns1_df) <- rnames
h5n1_ns1_df <- h5n1_ns1_df %>% mutate(Gene = "NS1", Subtype = "H5N1")

# Combine it all into one table
H5N1_all <- rbind(h5n1_pb2_df, h5n1_pb1_df, h5n1_pa_df, h5n1_ha_df, h5n1_np_df, h5n1_na_df, h5n1_m1_df, h5n1_ns1_df)

# Manipulate it to reflect reference species
long_H5N1 <- H5N1_all %>%
  pivot_longer(cols = starts_with("CAI"), names_to ="Reference Species", values_to = "CAI")

long_H5N1 <- long_H5N1 %>% mutate(Subtype = "H5N1")

################################################################################

#################################### H5N2 ######################################

################################################################################

# Import sequence files:
H5N2_pb2_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N2_genes\\PB2.fasta"
)
H5N2_pb1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N2_genes\\PB1.fasta"
)
H5N2_pa_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N2_genes\\polymerase PA.fasta"
)
H5N2_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N2_genes\\hemagglutinin.fasta"
)
H5N2_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N2_genes\\nucleocapsid.fasta"
)
H5N2_na_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N2_genes\\neuraminidase.fasta"
)
H5N2_m1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N2_genes\\matrix protein 1.fasta"
)
H5N2_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N2_genes\\nonstructural protein 1.fasta"
)

# Convert to codon tables
H5N2_pb2 <- codonTable(H5N2_pb2_dna)
H5N2_pb1 <- codonTable(H5N2_pb1_dna)
H5N2_pa <- codonTable(H5N2_pa_dna)
H5N2_ha <- codonTable(H5N2_ha_dna)
H5N2_np <- codonTable(H5N2_np_dna)
H5N2_na <- codonTable(H5N2_na_dna)
H5N2_m1 <- codonTable(H5N2_m1_dna)
H5N2_ns1 <- codonTable(H5N2_ns1_dna)

# Run the CAI analyses, convert to dataframes, manipulate df columns
h5n2_pb2_cai <- CAI(H5N2_pb2, subsets = list(chicken_ct, human_ct))
h5n2_pb2_df <- as.data.frame(h5n2_pb2_cai)
colnames(h5n2_pb2_df) <- rnames
h5n2_pb2_df <- h5n2_pb2_df %>% mutate(Gene = "PB2", Subtype = "H5N2")

h5n2_pb1_cai <- CAI(H5N2_pb1, subsets = list(chicken_ct, human_ct))
h5n2_pb1_df <- as.data.frame(h5n2_pb1_cai)
colnames(h5n2_pb1_df) <- rnames
h5n2_pb1_df <- h5n2_pb1_df %>% mutate(Gene = "PB1", Subtype = "H5N2")

h5n2_pa_cai <- CAI(H5N2_pa, subsets = list(chicken_ct, human_ct))
h5n2_pa_df <- as.data.frame(h5n2_pa_cai)
colnames(h5n2_pa_df) <- rnames
h5n2_pa_df <- h5n2_pa_df %>% mutate(Gene = "PA", Subtype = "H5N2")

h5n2_ha_cai <- CAI(H5N2_ha, subsets = list(chicken_ct, human_ct))
h5n2_ha_df <- as.data.frame(h5n2_ha_cai)
colnames(h5n2_ha_df) <- rnames
h5n2_ha_df <- h5n2_ha_df %>% mutate(Gene = "HA", Subtype = "H5N2")

h5n2_np_cai <- CAI(H5N2_np, subsets = list(chicken_ct, human_ct))
h5n2_np_df <- as.data.frame(h5n2_np_cai)
colnames(h5n2_np_df) <- rnames
h5n2_np_df <- h5n2_np_df %>% mutate(Gene = "NP", Subtype = "H5N2")

h5n2_na_cai <- CAI(H5N2_na, subsets = list(chicken_ct, human_ct))
h5n2_na_df <- as.data.frame(h5n2_na_cai)
colnames(h5n2_na_df) <- rnames
h5n2_na_df <- h5n2_na_df %>% mutate(Gene = "NA", Subtype = "H5N2")

h5n2_m1_cai <- CAI(H5N2_m1, subsets = list(chicken_ct, human_ct))
h5n2_m1_df <- as.data.frame(h5n2_m1_cai)
colnames(h5n2_m1_df) <- rnames
h5n2_m1_df <- h5n2_m1_df %>% mutate(Gene = "M1", Subtype = "H5N2")

h5n2_ns1_cai <- CAI(H5N2_ns1, subsets = list(chicken_ct, human_ct))
h5n2_ns1_df <- as.data.frame(h5n2_ns1_cai)
colnames(h5n2_ns1_df) <- rnames
h5n2_ns1_df <- h5n2_ns1_df %>% mutate(Gene = "NS1", Subtype = "H5N2")

# Combine it all into one table
H5N2_all <- rbind(h5n2_pb2_df, h5n2_pb1_df, h5n2_pa_df, h5n2_ha_df, h5n2_np_df, h5n2_na_df, h5n2_m1_df, h5n2_ns1_df)

# Manipulate it to reflect reference species
long_H5N2 <- H5N2_all %>%
  pivot_longer(cols = starts_with("CAI"), names_to ="Reference Species", values_to = "CAI")

long_H5N2 <- long_H5N2 %>% mutate(Subtype = "H5N2")

################################################################################

#################################### H5N3 ######################################

################################################################################

# Import sequence files:
H5N3_pb2_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N3_genes\\PB2.fasta"
)
H5N3_pb1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N3_genes\\PB1.fasta"
)
H5N3_pa_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N3_genes\\PA.fasta"
)
H5N3_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N3_genes\\HA.fasta"
)
H5N3_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N3_genes\\NP.fasta"
)
H5N3_na_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N3_genes\\NA.fasta"
)
H5N3_m1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N3_genes\\M1.fasta"
)
H5N3_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N3_genes\\NS1.fasta"
)

# Convert to codon tables
H5N3_pb2 <- codonTable(H5N3_pb2_dna)
H5N3_pb1 <- codonTable(H5N3_pb1_dna)
H5N3_pa <- codonTable(H5N3_pa_dna)
H5N3_ha <- codonTable(H5N3_ha_dna)
H5N3_np <- codonTable(H5N3_np_dna)
H5N3_na <- codonTable(H5N3_na_dna)
H5N3_m1 <- codonTable(H5N3_m1_dna)
H5N3_ns1 <- codonTable(H5N3_ns1_dna)

# Establish column names for the data frame
rnames <- c("CAI-Chicken", "CAI-Human")

# Run the CAI analyses, convert to dataframes, manipulate df columns
h5n3_pb2_cai <- CAI(H5N3_pb2, subsets = list(chicken_ct, human_ct))
h5n3_pb2_df <- as.data.frame(h5n3_pb2_cai)
colnames(h5n3_pb2_df) <- rnames
h5n3_pb2_df <- h5n3_pb2_df %>% mutate(Gene = "PB2", Subtype = "H5N3") #Add columns for the boxplot later

h5n3_pb1_cai <- CAI(H5N3_pb1, subsets = list(chicken_ct, human_ct))
h5n3_pb1_df <- as.data.frame(h5n3_pb1_cai)
colnames(h5n3_pb1_df) <- rnames
h5n3_pb1_df <- h5n3_pb1_df %>% mutate(Gene = "PB1", Subtype = "H5N3")

h5n3_pa_cai <- CAI(H5N3_pa, subsets = list(chicken_ct, human_ct))
h5n3_pa_df <- as.data.frame(h5n3_pa_cai)
colnames(h5n3_pa_df) <- rnames
h5n3_pa_df <- h5n3_pa_df %>% mutate(Gene = "PA", Subtype = "H5N3")

h5n3_ha_cai <- CAI(H5N3_ha, subsets = list(chicken_ct, human_ct))
h5n3_ha_df <- as.data.frame(h5n3_ha_cai)
colnames(h5n3_ha_df) <- rnames
h5n3_ha_df <- h5n3_ha_df %>% mutate(Gene = "HA", Subtype = "H5N3")

h5n3_np_cai <- CAI(H5N3_np, subsets = list(chicken_ct, human_ct))
h5n3_np_df <- as.data.frame(h5n3_np_cai)
colnames(h5n3_np_df) <- rnames
h5n3_np_df <- h5n3_np_df %>% mutate(Gene = "NP", Subtype = "H5N3")

h5n3_na_cai <- CAI(H5N3_na, subsets = list(chicken_ct, human_ct))
h5n3_na_df <- as.data.frame(h5n3_na_cai)
colnames(h5n3_na_df) <- rnames
h5n3_na_df <- h5n3_na_df %>% mutate(Gene = "NA", Subtype = "H5N3")

h5n3_m1_cai <- CAI(H5N3_m1, subsets = list(chicken_ct, human_ct))
h5n3_m1_df <- as.data.frame(h5n3_m1_cai)
colnames(h5n3_m1_df) <- rnames
h5n3_m1_df <- h5n3_m1_df %>% mutate(Gene = "M1", Subtype = "H5N3")

h5n3_ns1_cai <- CAI(H5N3_ns1, subsets = list(chicken_ct, human_ct))
h5n3_ns1_df <- as.data.frame(h5n3_ns1_cai)
colnames(h5n3_ns1_df) <- rnames
h5n3_ns1_df <- h5n3_ns1_df %>% mutate(Gene = "NS1", Subtype = "H5N3")

# Combine it all into one table
H5N3_all <- rbind(h5n3_pb2_df, h5n3_pb1_df, h5n3_pa_df, h5n3_ha_df, h5n3_np_df, h5n3_na_df, h5n3_m1_df, h5n3_ns1_df)

# Manipulate it to reflect reference species
long_H5N3 <- H5N3_all %>%
  pivot_longer(cols = starts_with("CAI"), names_to ="Reference Species", values_to = "CAI")

long_H5N3 <- long_H5N3 %>% mutate(Subtype = "H5N3")

################################################################################

################################## H5N4 ########################################

################################################################################

h5n4_pb2_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N4_genes\\PB2.fasta"
)
h5n4_pb2 <- codonTable(h5n4_pb2_dna)

h5n4_pb1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N4_genes\\PB1.fasta"
)
h5n4_pb1 <- codonTable(h5n4_pb1_dna)

h5n4_pa_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N4_genes\\PA.fasta"
)
h5n4_pa <- codonTable(h5n4_pa_dna)

h5n4_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N4_genes\\hemagglutinin.fasta"
)
h5n4_ha <- codonTable(h5n4_ha_dna)

h5n4_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N4_genes\\nucleocapsid.fasta"
)
h5n4_np <- codonTable(h5n4_np_dna)

h5n4_na_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N2_genes\\neuraminidase.fasta"
)
h5n4_na <- codonTable(h5n4_na_dna)

h5n4_mp1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N4_genes\\matrix protein 1.fasta"
)
h5n4_mp1 <- codonTable(h5n4_mp1_dna)

h5n4_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N4_genes\\nonstructural protein 1.fasta"
)
h5n4_ns1 <- codonTable(h5n4_ns1_dna)

# Run the CAI analyses, convert to dataframes, manipulate df columns
h5n4_pb2_cai <- CAI(h5n4_pb2, subsets = list(chicken_ct, human_ct))
h5n4_pb2_df <- as.data.frame(h5n4_pb2_cai)
colnames(h5n4_pb2_df) <- rnames
h5n4_pb2_df <- h5n4_pb2_df %>% mutate(Gene = "PB2", Subtype = "H5N4")

h5n4_pb1_cai <- CAI(h5n4_pb1, subsets = list(chicken_ct, human_ct))
h5n4_pb1_df <- as.data.frame(h5n4_pb1_cai)
colnames(h5n4_pb1_df) <- rnames
h5n4_pb1_df <- h5n4_pb1_df %>% mutate(Gene = "PB1", Subtype = "H5N4")

h5n4_pa_cai <- CAI(h5n4_pa, subsets = list(chicken_ct, human_ct))
h5n4_pa_df <- as.data.frame(h5n4_pa_cai)
colnames(h5n4_pa_df) <- rnames
h5n4_pa_df <- h5n4_pa_df %>% mutate(Gene = "PA", Subtype = "H5N4")

h5n4_ha_cai <- CAI(h5n4_ha, subsets = list(chicken_ct, human_ct))
h5n4_ha_df <- as.data.frame(h5n4_ha_cai)
colnames(h5n4_ha_df) <- rnames
h5n4_ha_df <- h5n4_ha_df %>% mutate(Gene = "HA", Subtype = "H5N4")

h5n4_np_cai <- CAI(h5n4_np, subsets = list(chicken_ct, human_ct))
h5n4_np_df <- as.data.frame(h5n4_np_cai)
colnames(h5n4_np_df) <- rnames
h5n4_np_df <- h5n4_np_df %>% mutate(Gene = "NP", Subtype = "H5N4")

h5n4_na_cai <- CAI(h5n4_na, subsets = list(chicken_ct, human_ct))
h5n4_na_df <- as.data.frame(h5n4_na_cai)
colnames(h5n4_na_df) <- rnames
h5n4_na_df <- h5n4_na_df %>% mutate(Gene = "NA", Subtype = "H5N4")

h5n4_m1_cai <- CAI(h5n4_mp1, subsets = list(chicken_ct, human_ct))
h5n4_m1_df <- as.data.frame(h5n4_m1_cai)
colnames(h5n4_m1_df) <- rnames
h5n4_m1_df <- h5n4_m1_df %>% mutate(Gene = "M1", Subtype = "H5N4")

h5n4_ns1_cai <- CAI(h5n4_ns1, subsets = list(chicken_ct, human_ct))
h5n4_ns1_df <- as.data.frame(h5n4_ns1_cai)
colnames(h5n4_ns1_df) <- rnames
h5n4_ns1_df <- h5n4_ns1_df %>% mutate(Gene = "NS1", Subtype = "H5N4")

# Combine it all into one table
H5N4_all <- rbind(h5n4_pb2_df, h5n4_pb1_df, h5n4_pa_df, h5n4_ha_df, h5n4_np_df, h5n4_na_df, h5n4_m1_df, h5n4_ns1_df)

# Manipulate it to reflect reference species
long_H5N4 <- H5N4_all %>%
  pivot_longer(cols = starts_with("CAI"), names_to ="Reference Species", values_to = "CAI")

long_H5N4 <- long_H5N4 %>% mutate(Subtype = "H5N4")

################################################################################

################################# H5N5 #########################################

################################################################################

h5n5_pb2_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N5_genes\\PB2.fasta"
)
h5n5_pb2 <- codonTable(h5n5_pb2_dna)

h5n5_pb1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N5_genes\\PB1.fasta"
)
h5n5_pb1 <- codonTable(h5n5_pb1_dna)

h5n5_pa_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N5_genes\\PA.fasta"
)
h5n5_pa <- codonTable(h5n5_pa_dna)

h5n5_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N5_genes\\hemagglutinin.fasta"
)
h5n5_ha <- codonTable(h5n5_ha_dna)

h5n5_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N5_genes\\nucleocapsid.fasta"
)
h5n5_np <- codonTable(h5n5_np_dna)

h5n5_na_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N5_genes\\neuraminidase.fasta"
)
h5n5_na <- codonTable(h5n5_na_dna)

h5n5_mp1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N5_genes\\matrix protein 1.fasta"
)
h5n5_mp1 <- codonTable(h5n5_mp1_dna)

h5n5_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N5_genes\\nonstructural protein 1.fasta"
)
h5n5_ns1 <- codonTable(h5n5_ns1_dna)

# Run the CAI analyses, convert to dataframes, manipulate df columns
h5n5_pb2_cai <- CAI(h5n5_pb2, subsets = list(chicken_ct, human_ct))
h5n5_pb2_df <- as.data.frame(h5n5_pb2_cai)
colnames(h5n5_pb2_df) <- rnames
h5n5_pb2_df <- h5n5_pb2_df %>% mutate(Gene = "PB2", Subtype = "H5N5")

h5n5_pb1_cai <- CAI(h5n5_pb1, subsets = list(chicken_ct, human_ct))
h5n5_pb1_df <- as.data.frame(h5n5_pb1_cai)
colnames(h5n5_pb1_df) <- rnames
h5n5_pb1_df <- h5n5_pb1_df %>% mutate(Gene = "PB1", Subtype = "H5N5")

h5n5_pa_cai <- CAI(h5n5_pa, subsets = list(chicken_ct, human_ct))
h5n5_pa_df <- as.data.frame(h5n5_pa_cai)
colnames(h5n5_pa_df) <- rnames
h5n5_pa_df <- h5n5_pa_df %>% mutate(Gene = "PA", Subtype = "H5N5")

h5n5_ha_cai <- CAI(h5n5_ha, subsets = list(chicken_ct, human_ct))
h5n5_ha_df <- as.data.frame(h5n5_ha_cai)
colnames(h5n5_ha_df) <- rnames
h5n5_ha_df <- h5n5_ha_df %>% mutate(Gene = "HA", Subtype = "H5N5")

h5n5_np_cai <- CAI(h5n5_np, subsets = list(chicken_ct, human_ct))
h5n5_np_df <- as.data.frame(h5n5_np_cai)
colnames(h5n5_np_df) <- rnames
h5n5_np_df <- h5n5_np_df %>% mutate(Gene = "NP", Subtype = "H5N5")

h5n5_na_cai <- CAI(h5n5_na, subsets = list(chicken_ct, human_ct))
h5n5_na_df <- as.data.frame(h5n5_na_cai)
colnames(h5n5_na_df) <- rnames
h5n5_na_df <- h5n5_na_df %>% mutate(Gene = "NA", Subtype = "H5N5")

h5n5_m1_cai <- CAI(h5n5_mp1, subsets = list(chicken_ct, human_ct))
h5n5_m1_df <- as.data.frame(h5n5_m1_cai)
colnames(h5n5_m1_df) <- rnames
h5n5_m1_df <- h5n5_m1_df %>% mutate(Gene = "M1", Subtype = "H5N5")

h5n5_ns1_cai <- CAI(h5n5_ns1, subsets = list(chicken_ct, human_ct))
h5n5_ns1_df <- as.data.frame(h5n5_ns1_cai)
colnames(h5n5_ns1_df) <- rnames
h5n5_ns1_df <- h5n5_ns1_df %>% mutate(Gene = "NS1", Subtype = "H5N5")

# Combine it all into one table
H5N5_all <- rbind(h5n5_pb2_df, h5n5_pb1_df, h5n5_pa_df, h5n5_ha_df, h5n5_np_df, h5n5_na_df, h5n5_m1_df, h5n5_ns1_df)

# Manipulate it to reflect reference species
long_H5N5 <- H5N5_all %>%
  pivot_longer(cols = starts_with("CAI"), names_to ="Reference Species", values_to = "CAI")

long_H5N5 <- long_H5N5 %>% mutate(Subtype = "H5N5")

##############################################################################

####################### All H5N6 sequences ###############################

##############################################################################

h5n6_pb2_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N6_genes\\PB2.fasta"
)
h5n6_pb2 <- codonTable(h5n6_pb2_dna)

h5n6_pb1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N6_genes\\PB1.fasta"
)
h5n6_pb1 <- codonTable(h5n6_pb1_dna)

h5n6_pa_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N6_genes\\PA.fasta"
)
h5n6_pa <- codonTable(h5n6_pa_dna)

h5n6_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N6_genes\\hemagglutinin.fasta"
)
h5n6_ha <- codonTable(h5n6_ha_dna)

h5n6_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N6_genes\\nucleocapsid.fasta"
)
h5n6_np <- codonTable(h5n6_np_dna)

h5n6_na_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N6_genes\\neuraminidase.fasta"
)
h5n6_na <- codonTable(h5n6_na_dna)

h5n6_mp1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N6_genes\\matrix protein 1.fasta"
)
h5n6_mp1 <- codonTable(h5n6_mp1_dna)

h5n6_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N6_genes\\nonstructural protein 1.fasta"
)
h5n6_ns1 <- codonTable(h5n6_ns1_dna)

# Run the CAI analyses, convert to dataframes, manipulate df columns
h5n6_pb2_cai <- CAI(h5n6_pb2, subsets = list(chicken_ct, human_ct))
h5n6_pb2_df <- as.data.frame(h5n6_pb2_cai)
colnames(h5n6_pb2_df) <- rnames
h5n6_pb2_df <- h5n6_pb2_df %>% mutate(Gene = "PB2", Subtype = "H5N6")

h5n6_pb1_cai <- CAI(h5n6_pb1, subsets = list(chicken_ct, human_ct))
h5n6_pb1_df <- as.data.frame(h5n6_pb1_cai)
colnames(h5n6_pb1_df) <- rnames
h5n6_pb1_df <- h5n6_pb1_df %>% mutate(Gene = "PB1", Subtype = "H5N6")

h5n6_pa_cai <- CAI(h5n6_pa, subsets = list(chicken_ct, human_ct))
h5n6_pa_df <- as.data.frame(h5n6_pa_cai)
colnames(h5n6_pa_df) <- rnames
h5n6_pa_df <- h5n6_pa_df %>% mutate(Gene = "PA", Subtype = "H5N6")

h5n6_ha_cai <- CAI(h5n6_ha, subsets = list(chicken_ct, human_ct))
h5n6_ha_df <- as.data.frame(h5n6_ha_cai)
colnames(h5n6_ha_df) <- rnames
h5n6_ha_df <- h5n6_ha_df %>% mutate(Gene = "HA", Subtype = "H5N6")

h5n6_np_cai <- CAI(h5n6_np, subsets = list(chicken_ct, human_ct))
h5n6_np_df <- as.data.frame(h5n6_np_cai)
colnames(h5n6_np_df) <- rnames
h5n6_np_df <- h5n6_np_df %>% mutate(Gene = "NP", Subtype = "H5N6")

h5n6_na_cai <- CAI(h5n6_na, subsets = list(chicken_ct, human_ct))
h5n6_na_df <- as.data.frame(h5n6_na_cai)
colnames(h5n6_na_df) <- rnames
h5n6_na_df <- h5n6_na_df %>% mutate(Gene = "NA", Subtype = "H5N6")

h5n6_m1_cai <- CAI(h5n6_mp1, subsets = list(chicken_ct, human_ct))
h5n6_m1_df <- as.data.frame(h5n6_m1_cai)
colnames(h5n6_m1_df) <- rnames
h5n6_m1_df <- h5n6_m1_df %>% mutate(Gene = "M1", Subtype = "H5N6")

h5n6_ns1_cai <- CAI(h5n6_ns1, subsets = list(chicken_ct, human_ct))
h5n6_ns1_df <- as.data.frame(h5n6_ns1_cai)
colnames(h5n6_ns1_df) <- rnames
h5n6_ns1_df <- h5n6_ns1_df %>% mutate(Gene = "NS1", Subtype = "H5N6")

# Combine it all into one table
H5N6_all <- rbind(h5n6_pb2_df, h5n6_pb1_df, h5n6_pa_df, h5n6_ha_df, h5n6_np_df, h5n6_na_df, h5n6_m1_df, h5n6_ns1_df)

# Manipulate it to reflect reference species
long_H5N6 <- H5N6_all %>%
  pivot_longer(cols = starts_with("CAI"), names_to ="Reference Species", values_to = "CAI")

long_H5N6 <- long_H5N6 %>% mutate(Subtype = "H5N6")

##############################################################################

################################## H5N7 #####################################

##############################################################################

h5n7_pb2_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N7_genes\\PB2.fasta"
)
h5n7_pb2 <- codonTable(h5n7_pb2_dna)

h5n7_pb1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N7_genes\\PB1.fasta"
)
h5n7_pb1 <- codonTable(h5n7_pb1_dna)

h5n7_pa_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N7_genes\\PA.fasta"
)
h5n7_pa <- codonTable(h5n7_pa_dna)

h5n7_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N7_genes\\hemagglutinin.fasta"
)
h5n7_ha <- codonTable(h5n7_ha_dna)

h5n7_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N7_genes\\nucleocapsid.fasta"
)
h5n7_np <- codonTable(h5n7_np_dna)

h5n7_na_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N7_genes\\neuraminidase.fasta"
)
h5n7_na <- codonTable(h5n7_na_dna)

h5n7_mp1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N7_genes\\matrix protein 1.fasta"
)
h5n7_mp1 <- codonTable(h5n7_mp1_dna)

h5n7_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N7_genes\\nonstructural protein 1.fasta"
)
h5n7_ns1 <- codonTable(h5n7_ns1_dna)

# Run the CAI analyses, convert to dataframes, manipulate df columns
h5n7_pb2_cai <- CAI(h5n7_pb2, subsets = list(chicken_ct, human_ct))
h5n7_pb2_df <- as.data.frame(h5n7_pb2_cai)
colnames(h5n7_pb2_df) <- rnames
h5n7_pb2_df <- h5n7_pb2_df %>% mutate(Gene = "PB2", Subtype = "H5N7")

h5n7_pb1_cai <- CAI(h5n7_pb1, subsets = list(chicken_ct, human_ct))
h5n7_pb1_df <- as.data.frame(h5n7_pb1_cai)
colnames(h5n7_pb1_df) <- rnames
h5n7_pb1_df <- h5n7_pb1_df %>% mutate(Gene = "PB1", Subtype = "H5N7")

h5n7_pa_cai <- CAI(h5n7_pa, subsets = list(chicken_ct, human_ct))
h5n7_pa_df <- as.data.frame(h5n7_pa_cai)
colnames(h5n7_pa_df) <- rnames
h5n7_pa_df <- h5n7_pa_df %>% mutate(Gene = "PA", Subtype = "H5N7")

h5n7_ha_cai <- CAI(h5n7_ha, subsets = list(chicken_ct, human_ct))
h5n7_ha_df <- as.data.frame(h5n7_ha_cai)
colnames(h5n7_ha_df) <- rnames
h5n7_ha_df <- h5n7_ha_df %>% mutate(Gene = "HA", Subtype = "H5N7")

h5n7_np_cai <- CAI(h5n7_np, subsets = list(chicken_ct, human_ct))
h5n7_np_df <- as.data.frame(h5n7_np_cai)
colnames(h5n7_np_df) <- rnames
h5n7_np_df <- h5n7_np_df %>% mutate(Gene = "NP", Subtype = "H5N7")

h5n7_na_cai <- CAI(h5n7_na, subsets = list(chicken_ct, human_ct))
h5n7_na_df <- as.data.frame(h5n7_na_cai)
colnames(h5n7_na_df) <- rnames
h5n7_na_df <- h5n7_na_df %>% mutate(Gene = "NA", Subtype = "H5N7")

h5n7_m1_cai <- CAI(h5n7_mp1, subsets = list(chicken_ct, human_ct))
h5n7_m1_df <- as.data.frame(h5n7_m1_cai)
colnames(h5n7_m1_df) <- rnames
h5n7_m1_df <- h5n7_m1_df %>% mutate(Gene = "M1", Subtype = "H5N7")

h5n7_ns1_cai <- CAI(h5n7_ns1, subsets = list(chicken_ct, human_ct))
h5n7_ns1_df <- as.data.frame(h5n7_ns1_cai)
colnames(h5n7_ns1_df) <- rnames
h5n7_ns1_df <- h5n7_ns1_df %>% mutate(Gene = "NS1", Subtype = "H5N7")

# Combine it all into one table
H5N7_all <- rbind(h5n7_pb2_df, h5n7_pb1_df, h5n7_pa_df, h5n7_ha_df, h5n7_np_df, h5n7_na_df, h5n7_m1_df, h5n7_ns1_df)

# Manipulate it to reflect reference species
long_H5N7 <- H5N7_all %>%
  pivot_longer(cols = starts_with("CAI"), names_to ="Reference Species", values_to = "CAI")

long_H5N7 <- long_H5N7 %>% mutate(Subtype = "H5N7")

##############################################################################

####################### All H5N8 sequences ###############################

##############################################################################

h5n8_pb2_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N8_genes\\PB2.fasta"
)
h5n8_pb2 <- codonTable(h5n8_pb2_dna)

h5n8_pb1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N8_genes\\PB1.fasta"
)
h5n8_pb1 <- codonTable(h5n8_pb1_dna)

h5n8_pa_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N8_genes\\PA.fasta"
)
h5n8_pa <- codonTable(h5n8_pa_dna)

h5n8_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N8_genes\\hemagglutinin.fasta"
)
h5n8_ha <- codonTable(h5n8_ha_dna)

h5n8_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N8_genes\\nucleocapsid.fasta"
)
h5n8_np <- codonTable(h5n8_np_dna)

h5n8_na_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N8_genes\\neuraminidase.fasta"
)
h5n8_na <- codonTable(h5n8_na_dna)

h5n8_mp1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N8_genes\\matrix protein 1.fasta"
)
h5n8_mp1 <- codonTable(h5n8_mp1_dna)

h5n8_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N8_genes\\nonstructural protein 1.fasta"
)
h5n8_ns1 <- codonTable(h5n8_ns1_dna)

# Run the CAI analyses, convert to dataframes, manipulate df columns
h5n8_pb2_cai <- CAI(h5n8_pb2, subsets = list(chicken_ct, human_ct))
h5n8_pb2_df <- as.data.frame(h5n8_pb2_cai)
colnames(h5n8_pb2_df) <- rnames
h5n8_pb2_df <- h5n8_pb2_df %>% mutate(Gene = "PB2", Subtype = "H5N8")

h5n8_pb1_cai <- CAI(h5n8_pb1, subsets = list(chicken_ct, human_ct))
h5n8_pb1_df <- as.data.frame(h5n8_pb1_cai)
colnames(h5n8_pb1_df) <- rnames
h5n8_pb1_df <- h5n8_pb1_df %>% mutate(Gene = "PB1", Subtype = "H5N8")

h5n8_pa_cai <- CAI(h5n8_pa, subsets = list(chicken_ct, human_ct))
h5n8_pa_df <- as.data.frame(h5n8_pa_cai)
colnames(h5n8_pa_df) <- rnames
h5n8_pa_df <- h5n8_pa_df %>% mutate(Gene = "PA", Subtype = "H5N8")

h5n8_ha_cai <- CAI(h5n8_ha, subsets = list(chicken_ct, human_ct))
h5n8_ha_df <- as.data.frame(h5n8_ha_cai)
colnames(h5n8_ha_df) <- rnames
h5n8_ha_df <- h5n8_ha_df %>% mutate(Gene = "HA", Subtype = "H5N8")

h5n8_np_cai <- CAI(h5n8_np, subsets = list(chicken_ct, human_ct))
h5n8_np_df <- as.data.frame(h5n8_np_cai)
colnames(h5n8_np_df) <- rnames
h5n8_np_df <- h5n8_np_df %>% mutate(Gene = "NP", Subtype = "H5N8")

h5n8_na_cai <- CAI(h5n8_na, subsets = list(chicken_ct, human_ct))
h5n8_na_df <- as.data.frame(h5n8_na_cai)
colnames(h5n8_na_df) <- rnames
h5n8_na_df <- h5n8_na_df %>% mutate(Gene = "NA", Subtype = "H5N8")

h5n8_m1_cai <- CAI(h5n8_mp1, subsets = list(chicken_ct, human_ct))
h5n8_m1_df <- as.data.frame(h5n8_m1_cai)
colnames(h5n8_m1_df) <- rnames
h5n8_m1_df <- h5n8_m1_df %>% mutate(Gene = "M1", Subtype = "H5N8")

h5n8_ns1_cai <- CAI(h5n8_ns1, subsets = list(chicken_ct, human_ct))
h5n8_ns1_df <- as.data.frame(h5n8_ns1_cai)
colnames(h5n8_ns1_df) <- rnames
h5n8_ns1_df <- h5n8_ns1_df %>% mutate(Gene = "NS1", Subtype = "H5N8")

# Combine it all into one table
H5N8_all <- rbind(h5n8_pb2_df, h5n8_pb1_df, h5n8_pa_df, h5n8_ha_df, h5n8_np_df, h5n8_na_df, h5n8_m1_df, h5n8_ns1_df)

# Manipulate it to reflect reference species
long_H5N8 <- H5N8_all %>%
  pivot_longer(cols = starts_with("CAI"), names_to ="Reference Species", values_to = "CAI")

long_H5N8 <- long_H5N8 %>% mutate(Subtype = "H5N8")

##############################################################################

####################### All H5N9 sequences ###############################

##############################################################################

h5n9_pb2_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N9_genes\\PB2.fasta"
)
h5n9_pb2 <- codonTable(h5n9_pb2_dna)

h5n9_pb1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N9_genes\\PB1.fasta"
)
h5n9_pb1 <- codonTable(h5n9_pb1_dna)

h5n9_pa_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N9_genes\\PA.fasta"
)
h5n9_pa <- codonTable(h5n9_pa_dna)

h5n9_ha_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N9_genes\\HA.fasta"
)
h5n9_ha <- codonTable(h5n9_ha_dna)

h5n9_np_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N9_genes\\NP.fasta"
)
h5n9_np <- codonTable(h5n9_np_dna)

h5n9_na_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N9_genes\\NA.fasta"
)
h5n9_na <- codonTable(h5n9_na_dna)

h5n9_mp1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N9_genes\\M1.fasta"
)
h5n9_mp1 <- codonTable(h5n9_mp1_dna)

h5n9_ns1_dna <- readSet(
  file = "C:\\Users\\reedm\\OneDrive\\Documents\\BIOT 670\\Project\\Sequence Data\\H5N9_genes\\NS1.fasta"
)
h5n9_ns1 <- codonTable(h5n9_ns1_dna)

# Run the CAI analyses, convert to dataframes, manipulate df columns
h5n9_pb2_cai <- CAI(h5n9_pb2, subsets = list(chicken_ct, human_ct))
h5n9_pb2_df <- as.data.frame(h5n9_pb2_cai)
colnames(h5n9_pb2_df) <- rnames
h5n9_pb2_df <- h5n9_pb2_df %>% mutate(Gene = "PB2", Subtype = "H5N9")

h5n9_pb1_cai <- CAI(h5n9_pb1, subsets = list(chicken_ct, human_ct))
h5n9_pb1_df <- as.data.frame(h5n9_pb1_cai)
colnames(h5n9_pb1_df) <- rnames
h5n9_pb1_df <- h5n9_pb1_df %>% mutate(Gene = "PB1", Subtype = "H5N9")

h5n9_pa_cai <- CAI(h5n9_pa, subsets = list(chicken_ct, human_ct))
h5n9_pa_df <- as.data.frame(h5n9_pa_cai)
colnames(h5n9_pa_df) <- rnames
h5n9_pa_df <- h5n9_pa_df %>% mutate(Gene = "PA", Subtype = "H5N9")

h5n9_ha_cai <- CAI(h5n9_ha, subsets = list(chicken_ct, human_ct))
h5n9_ha_df <- as.data.frame(h5n9_ha_cai)
colnames(h5n9_ha_df) <- rnames
h5n9_ha_df <- h5n9_ha_df %>% mutate(Gene = "HA", Subtype = "H5N9")

h5n9_np_cai <- CAI(h5n9_np, subsets = list(chicken_ct, human_ct))
h5n9_np_df <- as.data.frame(h5n9_np_cai)
colnames(h5n9_np_df) <- rnames
h5n9_np_df <- h5n5_np_df %>% mutate(Gene = "NP", Subtype = "H5N9")

h5n9_na_cai <- CAI(h5n9_na, subsets = list(chicken_ct, human_ct))
h5n9_na_df <- as.data.frame(h5n9_na_cai)
colnames(h5n9_na_df) <- rnames
h5n9_na_df <- h5n9_na_df %>% mutate(Gene = "NA", Subtype = "H5N9")

h5n9_m1_cai <- CAI(h5n9_mp1, subsets = list(chicken_ct, human_ct))
h5n9_m1_df <- as.data.frame(h5n9_m1_cai)
colnames(h5n9_m1_df) <- rnames
h5n9_m1_df <- h5n9_m1_df %>% mutate(Gene = "M1", Subtype = "H5N9")

h5n9_ns1_cai <- CAI(h5n9_ns1, subsets = list(chicken_ct, human_ct))
h5n9_ns1_df <- as.data.frame(h5n9_ns1_cai)
colnames(h5n9_ns1_df) <- rnames
h5n9_ns1_df <- h5n9_ns1_df %>% mutate(Gene = "NS1", Subtype = "H5N9")

# Combine it all into one table
H5N9_all <- rbind(h5n9_pb2_df, h5n9_pb1_df, h5n9_pa_df, h5n9_ha_df, h5n9_np_df, h5n9_na_df, h5n9_m1_df, h5n9_ns1_df)

# Manipulate it to reflect reference species
long_H5N9 <- H5N9_all %>%
  pivot_longer(cols = starts_with("CAI"), names_to ="Reference Species", values_to = "CAI")

long_H5N9 <- long_H5N9 %>% mutate(Subtype = "H5N9")

# Now merge all subtypes into one data frame
all_H5 <- bind_rows(long_H5N1, long_H5N2, long_H5N3, long_H5N4, long_H5N5, long_H5N6, long_H5N7, long_H5N8, long_H5N9)

# Plot it all
# View by Gene:
ggplot(all_H5, aes(x = Subtype, y = CAI, fill = `Reference Species`)) +
         geom_boxplot() +
         theme_classic() +
         labs(title = "H5 CAI Comparison by Gene", x = "Gene", y = "CAI") +
  facet_wrap(~ Gene) +
  geom_hline(yintercept = mean(all_H5$CAI), linetype = 2)


# By subtype
ggplot(all_H5, aes(x = Gene, y = CAI, fill = `Reference Species`)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "H5 CAI Comparison by Subtype", x = "Gene", y = "CAI") +
  facet_wrap(~ Subtype) +
  geom_hline(yintercept = mean(all_H5$CAI), linetype = 2)
  

# Run an anova to see if there's significant variation between genes by subtype
# ANOVAs for Chicken CAIs
# PB2
pb2_aov_df <- all_H5[all_H5$Gene == "PB2" & all_H5$`Reference Species` == "CAI-Chicken",]
pb2_aov <- aov(CAI ~ Subtype, pb2_aov_df)
summary(pb2_aov)

# PB1
pb1_aov_df <- all_H5[all_H5$Gene == "PB1" & all_H5$`Reference Species` == "CAI-Chicken",]
pb1_aov <- aov(CAI ~ Subtype, pb1_aov_df)
summary(pb1_aov)

# PA
pa_aov_df <- all_H5[all_H5$Gene == "PA" & all_H5$`Reference Species` == "CAI-Chicken",]
pa_aov <- aov(CAI ~ Subtype, pa_aov_df)
summary(pa_aov)

# HA
ha_aov_df <- all_H5[all_H5$Gene == "HA" & all_H5$`Reference Species` == "CAI-Chicken",]
ha_aov <- aov(CAI ~ Subtype, ha_aov_df)
summary(ha_aov)

# NP
np_aov_df <- all_H5[all_H5$Gene == "NP" & all_H5$`Reference Species` == "CAI-Chicken",]
np_aov <- aov(CAI ~ Subtype, np_aov_df)
summary(np_aov)

# NA
na_aov_df <- all_H5[all_H5$Gene == "NA" & all_H5$`Reference Species` == "CAI-Chicken",]
na_aov <- aov(CAI ~ Subtype, na_aov_df)
summary(na_aov)

# M1
m1_aov_df <- all_H5[all_H5$Gene == "M1" & all_H5$`Reference Species` == "CAI-Chicken",]
m1_aov <- aov(CAI ~ Subtype, m1_aov_df)
summary(m1_aov)

# NS1
ns1_aov_df <- all_H5[all_H5$Gene == "NS1" & all_H5$`Reference Species` == "CAI-Chicken",]
ns1_aov <- aov(CAI ~ Subtype, ns1_aov_df)
summary(ns1_aov)

# ANOVAs for human CAIs
# PB2
pb2_aovh_df <- all_H5[all_H5$Gene == "PB2" & all_H5$`Reference Species` == "CAI-Human",]
pb2_aovh <- aov(CAI ~ Subtype, pb2_aovh_df)
summary(pb2_aovh)

# PB1
pb1_aovh_df <- all_H5[all_H5$Gene == "PB1" & all_H5$`Reference Species` == "CAI-Human",]
pb1_aovh <- aov(CAI ~ Subtype, pb1_aovh_df)
summary(pb1_aovh)

# PA
pa_aovh_df <- all_H5[all_H5$Gene == "PA" & all_H5$`Reference Species` == "CAI-Human",]
pa_aovh <- aov(CAI ~ Subtype, pa_aovh_df)
summary(pa_aovh)

# HA
ha_aovh_df <- all_H5[all_H5$Gene == "HA" & all_H5$`Reference Species` == "CAI-Human",]
ha_aovh <- aov(CAI ~ Subtype, ha_aovh_df)
summary(ha_aovh)

# NP
np_aovh_df <- all_H5[all_H5$Gene == "NP" & all_H5$`Reference Species` == "CAI-Human",]
np_aovh <- aov(CAI ~ Subtype, np_aovh_df)
summary(np_aovh)

# NA
na_aovh_df <- all_H5[all_H5$Gene == "NA" & all_H5$`Reference Species` == "CAI-Human",]
na_aovh <- aov(CAI ~ Subtype, na_aovh_df)
summary(na_aovh)

# M1
m1_aovh_df <- all_H5[all_H5$Gene == "M1" & all_H5$`Reference Species` == "CAI-Human",]
m1_aovh <- aov(CAI ~ Subtype, m1_aovh_df)
summary(m1_aovh)

# NS1
ns1_aovh_df <- all_H5[all_H5$Gene == "NS1" & all_H5$`Reference Species` == "CAI-Human",]
ns1_aovh <- aov(CAI ~ Subtype, ns1_aovh_df)
summary(ns1_aovh)

