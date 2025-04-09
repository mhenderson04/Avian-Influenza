library(ggplot2)
library(tidyverse)
library(coRdon)

# Import sequence data

# H5N1
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

# H5N2
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

# H5N3
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

# H5N4
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

# H5N5
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

# H5N6
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

# H5N7
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

# H5N8
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

# H5N9
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

################################################################################

################################# ENc analysis #################################

################################################################################

#PB2 genes
h5n1_pb2_enc <- ENC(H5N1_pb2)
h5n1_pb2enc_df <- data.frame(Subtype = "H5N1",
                             ENC = h5n1_pb2_enc)
h5n2_pb2_enc <- ENC(H5N2_pb2)
h5n2_pb2enc_df <- data.frame(Subtype = "H5N2",
                             ENC = h5n2_pb2_enc)
h5n3_pb2_enc <- ENC(H5N3_pb2)
h5n3_pb2enc_df <- data.frame(Subtype = "H5N3",
                             ENC = h5n3_pb2_enc)
h5n4_pb2_enc <- ENC(h5n4_pb2)
h5n4_pb2enc_df <- data.frame(Subtype = "H5N4",
                             ENC = h5n4_pb2_enc)
h5n5_pb2_enc <- ENC(h5n5_pb2)
h5n5_pb2enc_df <- data.frame(Subtype = "H5N5",
                             ENC = h5n5_pb2_enc)
h5n6_pb2_enc <- ENC(h5n6_pb2)
h5n6_pb2enc_df <- data.frame(Subtype = "H5N6",
                             ENC = h5n6_pb2_enc)
h5n7_pb2_enc <- ENC(h5n7_pb2)
h5n7_pb2enc_df <- data.frame(Subtype = "H5N7",
                             ENC = h5n7_pb2_enc)
h5n8_pb2_enc <- ENC(h5n8_pb2)
h5n8_pb2enc_df <- data.frame(Subtype = "H5N8",
                             ENC = h5n8_pb2_enc)
h5n9_pb2_enc <- ENC(h5n9_pb2)
h5n9_pb2enc_df <- data.frame(Subtype = "H5N9",
                             ENC = h5n9_pb2_enc)

pb2_enc_df <- rbind(h5n1_pb2enc_df, h5n2_pb2enc_df, h5n3_pb2enc_df, h5n4_pb2enc_df, h5n5_pb2enc_df,
                    h5n6_pb2enc_df, h5n7_pb2enc_df, h5n8_pb2enc_df, h5n9_pb2enc_df)
pb2_enc_df <- pb2_enc_df %>% mutate(Gene = "PB2")

#PB1 genes
h5n1_pb1_enc <- ENC(H5N1_pb1)
h5n1_pb1enc_df <- data.frame(Subtype = "H5N1",
                             ENC = h5n1_pb1_enc)
h5n2_pb1_enc <- ENC(H5N2_pb1)
h5n2_pb1enc_df <- data.frame(Subtype = "H5N2",
                             ENC = h5n2_pb1_enc)
h5n3_pb1_enc <- ENC(H5N3_pb1)
h5n3_pb1enc_df <- data.frame(Subtype = "H5N3",
                             ENC = h5n3_pb1_enc)
h5n4_pb1_enc <- ENC(h5n4_pb1)
h5n4_pb1enc_df <- data.frame(Subtype = "H5N4",
                             ENC = h5n4_pb1_enc)
h5n5_pb1_enc <- ENC(h5n5_pb1)
h5n5_pb1enc_df <- data.frame(Subtype = "H5N5",
                             ENC = h5n5_pb1_enc)
h5n6_pb1_enc <- ENC(h5n6_pb1)
h5n6_pb1enc_df <- data.frame(Subtype = "H5N6",
                             ENC = h5n6_pb1_enc)
h5n7_pb1_enc <- ENC(h5n7_pb1)
h5n7_pb1enc_df <- data.frame(Subtype = "H5N7",
                             ENC = h5n7_pb1_enc)
h5n8_pb1_enc <- ENC(h5n8_pb1)
h5n8_pb1enc_df <- data.frame(Subtype = "H5N8",
                             ENC = h5n8_pb1_enc)
h5n9_pb1_enc <- ENC(h5n9_pb1)
h5n9_pb1enc_df <- data.frame(Subtype = "H5N9",
                             ENC = h5n9_pb1_enc)

pb1_enc_df <- rbind(h5n1_pb1enc_df, h5n2_pb1enc_df, h5n3_pb1enc_df, h5n4_pb1enc_df, h5n5_pb1enc_df,
                    h5n6_pb1enc_df, h5n7_pb1enc_df, h5n8_pb1enc_df, h5n9_pb1enc_df)
pb1_enc_df <- pb1_enc_df %>% mutate(Gene = "PB1")

#PA genes
h5n1_pa_enc <- ENC(H5N1_pa)
h5n1_paenc_df <- data.frame(Subtype = "H5N1",
                            ENC = h5n1_pa_enc)
h5n2_pa_enc <- ENC(H5N2_pa)
h5n2_paenc_df <- data.frame(Subtype = "H5N2",
                            ENC = h5n2_pa_enc)
h5n3_pa_enc <- ENC(H5N3_pa)
h5n3_paenc_df <- data.frame(Subtype = "H5N3",
                            ENC = h5n3_pa_enc)
h5n4_pa_enc <- ENC(h5n4_pa)
h5n4_paenc_df <- data.frame(Subtype = "H5N4",
                            ENC = h5n4_pa_enc)
h5n5_pa_enc <- ENC(h5n5_pa)
h5n5_paenc_df <- data.frame(Subtype = "H5N5",
                            ENC = h5n5_pa_enc)
h5n6_pa_enc <- ENC(h5n6_pa)
h5n6_paenc_df <- data.frame(Subtype = "H5N6",
                            ENC = h5n6_pa_enc)
h5n7_pa_enc <- ENC(h5n7_pa)
h5n7_paenc_df <- data.frame(Subtype = "H5N7",
                            ENC = h5n7_pa_enc)
h5n8_pa_enc <- ENC(h5n8_pa)
h5n8_paenc_df <- data.frame(Subtype = "H5N8",
                            ENC = h5n8_pa_enc)
h5n9_pa_enc <- ENC(h5n9_pa)
h5n9_paenc_df <- data.frame(Subtype = "H5N9",
                            ENC = h5n9_pa_enc)

pa_enc_df <- rbind(h5n1_paenc_df, h5n2_paenc_df, h5n3_paenc_df, h5n4_paenc_df, h5n5_paenc_df,
                   h5n6_paenc_df, h5n7_paenc_df, h5n8_paenc_df, h5n9_paenc_df)
pa_enc_df <- pa_enc_df %>% mutate(Gene = "PA")

# HA genes
h5n1_ha_enc <- ENC(H5N1_ha)
h5n1_haenc_df <- data.frame(Subtype = "H5N1",
                            ENC = h5n1_ha_enc)
h5n2_ha_enc <- ENC(H5N2_ha)
h5n2_haenc_df <- data.frame(Subtype = "H5N2",
                            ENC = h5n2_ha_enc)
h5n3_ha_enc <- ENC(H5N3_ha)
h5n3_haenc_df <- data.frame(Subtype = "H5N3",
                            ENC = h5n3_ha_enc)
h5n4_ha_enc <- ENC(h5n4_ha)
h5n4_haenc_df <- data.frame(Subtype = "H5N4",
                            ENC = h5n4_ha_enc)
h5n5_ha_enc <- ENC(h5n5_ha)
h5n5_haenc_df <- data.frame(Subtype = "H5N5",
                            ENC = h5n5_ha_enc)
h5n6_ha_enc <- ENC(h5n6_ha)
h5n6_haenc_df <- data.frame(Subtype = "H5N6",
                            ENC = h5n6_ha_enc)
h5n7_ha_enc <- ENC(h5n7_ha)
h5n7_haenc_df <- data.frame(Subtype = "H5N7",
                            ENC = h5n7_ha_enc)
h5n8_ha_enc <- ENC(h5n8_ha)
h5n8_haenc_df <- data.frame(Subtype = "H5N8",
                            ENC = h5n8_ha_enc)
h5n9_ha_enc <- ENC(h5n9_ha)
h5n9_haenc_df <- data.frame(Subtype = "H5N9",
                            ENC = h5n9_ha_enc)

ha_enc_df <- rbind(h5n1_haenc_df, h5n2_haenc_df, h5n3_haenc_df, h5n4_haenc_df, h5n5_haenc_df,
                   h5n6_haenc_df, h5n7_haenc_df, h5n8_haenc_df, h5n9_haenc_df)
ha_enc_df <- ha_enc_df %>% mutate(Gene = "HA")

# NP genes
h5n1_np_enc <- ENC(H5N1_np)
h5n1_npenc_df <- data.frame(Subtype = "H5N1",
                            ENC = h5n1_np_enc)
h5n2_np_enc <- ENC(H5N2_np)
h5n2_npenc_df <- data.frame(Subtype = "H5N2",
                            ENC = h5n2_np_enc)
h5n3_np_enc <- ENC(H5N3_np)
h5n3_npenc_df <- data.frame(Subtype = "H5N3",
                            ENC = h5n3_np_enc)
h5n4_np_enc <- ENC(h5n4_np)
h5n4_npenc_df <- data.frame(Subtype = "H5N4",
                            ENC = h5n4_np_enc)
h5n5_np_enc <- ENC(h5n5_np)
h5n5_npenc_df <- data.frame(Subtype = "H5N5",
                            ENC = h5n5_np_enc)
h5n6_np_enc <- ENC(h5n6_np)
h5n6_npenc_df <- data.frame(Subtype = "H5N6",
                            ENC = h5n6_np_enc)
h5n7_np_enc <- ENC(h5n7_np)
h5n7_npenc_df <- data.frame(Subtype = "H5N7",
                            ENC = h5n7_np_enc)
h5n8_np_enc <- ENC(h5n8_np)
h5n8_npenc_df <- data.frame(Subtype = "H5N8",
                            ENC = h5n8_np_enc)
h5n9_np_enc <- ENC(h5n9_np)
h5n9_npenc_df <- data.frame(Subtype = "H5N9",
                            ENC = h5n9_np_enc)

np_enc_df <- rbind(h5n1_npenc_df, h5n2_npenc_df, h5n3_npenc_df, h5n4_npenc_df, h5n5_npenc_df,
                   h5n6_npenc_df, h5n7_npenc_df, h5n8_npenc_df, h5n9_npenc_df)
np_enc_df <- np_enc_df %>% mutate(Gene = "NP")

# NA genes
h5n1_na_enc <- ENC(H5N1_na)
h5n1_naenc_df <- data.frame(Subtype = "H5N1",
                            ENC = h5n1_na_enc)
h5n2_na_enc <- ENC(H5N2_na)
h5n2_naenc_df <- data.frame(Subtype = "H5N2",
                            ENC = h5n2_na_enc)
h5n3_na_enc <- ENC(H5N3_na)
h5n3_naenc_df <- data.frame(Subtype = "H5N3",
                            ENC = h5n3_na_enc)
h5n4_na_enc <- ENC(h5n4_na)
h5n4_naenc_df <- data.frame(Subtype = "H5N4",
                            ENC = h5n4_na_enc)
h5n5_na_enc <- ENC(h5n5_na)
h5n5_naenc_df <- data.frame(Subtype = "H5N5",
                            ENC = h5n5_na_enc)
h5n6_na_enc <- ENC(h5n6_na)
h5n6_naenc_df <- data.frame(Subtype = "H5N6",
                            ENC = h5n6_na_enc)
h5n7_na_enc <- ENC(h5n7_na)
h5n7_naenc_df <- data.frame(Subtype = "H5N7",
                            ENC = h5n7_na_enc)
h5n8_na_enc <- ENC(h5n8_na)
h5n8_naenc_df <- data.frame(Subtype = "H5N8",
                            ENC = h5n8_na_enc)
h5n9_na_enc <- ENC(h5n9_na)
h5n9_naenc_df <- data.frame(Subtype = "H5N9",
                            ENC = h5n9_na_enc)

na_enc_df <- rbind(h5n1_naenc_df, h5n2_naenc_df, h5n3_naenc_df, h5n4_naenc_df, h5n5_naenc_df,
                   h5n6_naenc_df, h5n7_naenc_df, h5n8_naenc_df, h5n9_naenc_df)
na_enc_df <- na_enc_df %>% mutate(Gene = "NA")

# M1 genes
h5n1_m1_enc <- ENC(H5N1_m1)
h5n1_m1enc_df <- data.frame(Subtype = "H5N1",
                            ENC = h5n1_m1_enc)
h5n2_m1_enc <- ENC(H5N2_m1)
h5n2_m1enc_df <- data.frame(Subtype = "H5N2",
                            ENC = h5n2_m1_enc)
h5n3_m1_enc <- ENC(H5N3_m1)
h5n3_m1enc_df <- data.frame(Subtype = "H5N3",
                            ENC = h5n3_m1_enc)
h5n4_m1_enc <- ENC(h5n4_mp1)
h5n4_m1enc_df <- data.frame(Subtype = "H5N4",
                            ENC = h5n4_m1_enc)
h5n5_m1_enc <- ENC(h5n5_mp1)
h5n5_m1enc_df <- data.frame(Subtype = "H5N5",
                            ENC = h5n5_m1_enc)
h5n6_m1_enc <- ENC(h5n6_mp1)
h5n6_m1enc_df <- data.frame(Subtype = "H5N6",
                            ENC = h5n6_m1_enc)
h5n7_m1_enc <- ENC(h5n7_mp1)
h5n7_m1enc_df <- data.frame(Subtype = "H5N7",
                            ENC = h5n7_m1_enc)
h5n8_m1_enc <- ENC(h5n8_mp1)
h5n8_m1enc_df <- data.frame(Subtype = "H5N8",
                            ENC = h5n8_m1_enc)
h5n9_m1_enc <- ENC(h5n9_mp1)
h5n9_m1enc_df <- data.frame(Subtype = "H5N9",
                            ENC = h5n9_m1_enc)

m1_enc_df <- rbind(h5n1_m1enc_df, h5n2_m1enc_df, h5n3_m1enc_df, h5n4_m1enc_df, h5n5_m1enc_df,
                   h5n6_m1enc_df, h5n7_m1enc_df, h5n8_m1enc_df, h5n9_m1enc_df)
m1_enc_df <- m1_enc_df %>% mutate(Gene = "M1")

# NS1 genes
h5n1_ns1_enc <- ENC(H5N1_ns1)
h5n1_ns1enc_df <- data.frame(Subtype = "H5N1",
                             ENC = h5n1_ns1_enc)
h5n2_ns1_enc <- ENC(H5N2_ns1)
h5n2_ns1enc_df <- data.frame(Subtype = "H5N2",
                             ENC = h5n2_ns1_enc)
h5n3_ns1_enc <- ENC(H5N3_ns1)
h5n3_ns1enc_df <- data.frame(Subtype = "H5N3",
                             ENC = h5n3_ns1_enc)
h5n4_ns1_enc <- ENC(h5n4_ns1)
h5n4_ns1enc_df <- data.frame(Subtype = "H5N4",
                             ENC = h5n4_ns1_enc)
h5n5_ns1_enc <- ENC(h5n5_ns1)
h5n5_ns1enc_df <- data.frame(Subtype = "H5N5",
                             ENC = h5n5_ns1_enc)
h5n6_ns1_enc <- ENC(h5n6_ns1)
h5n6_ns1enc_df <- data.frame(Subtype = "H5N6",
                             ENC = h5n6_ns1_enc)
h5n7_ns1_enc <- ENC(h5n7_ns1)
h5n7_ns1enc_df <- data.frame(Subtype = "H5N7",
                             ENC = h5n7_ns1_enc)
h5n8_ns1_enc <- ENC(h5n8_ns1)
h5n8_ns1enc_df <- data.frame(Subtype = "H5N8",
                             ENC = h5n8_ns1_enc)
h5n9_ns1_enc <- ENC(h5n9_ns1)
h5n9_ns1enc_df <- data.frame(Subtype = "H5N9",
                             ENC = h5n9_ns1_enc)

ns1_enc_df <- rbind(h5n1_ns1enc_df, h5n2_ns1enc_df, h5n3_ns1enc_df, h5n4_ns1enc_df, h5n5_ns1enc_df,
                    h5n6_ns1enc_df, h5n7_ns1enc_df, h5n8_ns1enc_df, h5n9_ns1enc_df)
ns1_enc_df <- ns1_enc_df %>% mutate(Gene = "NS1")

all_enc <- rbind(pb2_enc_df, pb1_enc_df, pa_enc_df, ha_enc_df, np_enc_df, na_enc_df, m1_enc_df, ns1_enc_df)
summary(all_enc)

all_enc %>%
  ggplot(aes(ENC, fill = Gene)) +
  geom_density(alpha = 0.5)