# Load the required libraries
library(bio3d)
library(seqinr)

# 1. Retrieve protein sequences and structures
insulin.seq <- read.fasta("https://www.uniprot.org/uniprot/P01308.fasta")
insulin.pdb <- read.pdb("https://files.rcsb.org/view/9INS.pdb")

keratin.seq <- read.fasta("https://www.uniprot.org/uniprot/P04264.fasta")
keratin.pdb <- read.pdb("https://files.rcsb.org/view/3ASW.pdb")

# 2. Find amino acid composition and plot
calculate_aa_composition <- function(seq) {
  aa_count <- table(seq)
  aa_names <- names(aa_count)
  aa_freq <- aa_count / sum(aa_count)
  return(data.frame(Amino_Acid = aa_names, Frequency = aa_freq))
}

insulin_comp <- calculate_aa_composition(insulin.seq[[1]])
keratin_comp <- calculate_aa_composition(keratin.seq[[1]])

# Convert Frequency.Freq column to numeric
insulin_comp$Frequency.Freq <- as.numeric(as.character(insulin_comp$Frequency.Freq))
keratin_comp$Frequency.Freq <- as.numeric(as.character(keratin_comp$Frequency.Freq))

# Plot the barplots again
par(mfrow = c(1, 2))
barplot(sort(insulin_comp$Frequency.Freq), main = "Insulin Amino Acid Composition", xlab = "Amino Acid", ylab = "Count")
barplot(sort(keratin_comp$Frequency.Freq), main = "Keratin Amino Acid Composition", xlab = "Amino Acid", ylab = "Count")

# 3. Use BLAST to find 10 most similar structures
insulin.blast <- blast.pdb(insulin.pdb)
keratin.blast <- blast.pdb(keratin.pdb)

hits <- plot.blast(insulin.blast, top = 10)
hits <- plot.blast(keratin.blast, top = 10)

preprocess_structure <- function(pdb_file) {
  pdb <- read.pdb(pdb_file)
  
  # Remove non-protein atoms and keep only the protein chain (e.g., chain A)
  protein_chain <- "A"  # Adjust this if the protein chain is different
  trimmed_pdb <- trim.pdb(pdb, chain = protein_chain)
  
  return(trimmed_pdb)
}

insulin_trimmed <- preprocess_structure("9ins")
keratin_trimmed <- preprocess_structure("3asw")

# Perform normal mode analysis (NMA) on the preprocessed structures
insulin.modes <- nma(insulin_trimmed)
keratin.modes <- nma(keratin_trimmed)

# Plot normal mode analysis
plot.nma(insulin.modes)
plot.nma(keratin.modes)

library(rgl)
# Extract coordinates from the PDB structures
insulin_coords <- trim.pdb(insulin.pdb, inds = atom.select(insulin.pdb, "calpha"))$atom
keratin_coords <- trim.pdb(keratin.pdb, inds = atom.select(keratin.pdb, "calpha"))$atom

# Check if coordinate data is valid
print(insulin_coords)
print(keratin_coords)

open3d()

plot3d(insulin_coords$x, insulin_coords$y, insulin_coords$z, col = "blue", type = "s", size = 1)
text3d(insulin_coords$centroid[1], insulin_coords$centroid[2], insulin_coords$centroid[3], text = "Insulin", adj = c(-0.5, 0))
title3d("3D Visualization of Insulin and Amino Acids")
axes3d(labels = c("X-axis", "Y-axis", "Z-axis"))
rgl.viewpoint(theta = 120, phi = 20, fov = 60, zoom = 0.8)

open3d()

plot3d(keratin_coords$x, keratin_coords$y, keratin_coords$z, col = "red", type = "s", size = 1)
text3d(keratin_coords$centroid[1], keratin_coords$centroid[2], keratin_coords$centroid[3], text = "Keratin", adj = c(-0.5, 0))
title3d("3D Visualization of Keratin and Amino Acids")
axes3d(labels = c("X-axis", "Y-axis", "Z-axis"))
rgl.viewpoint(theta = 120, phi = 20, fov = 60, zoom = 0.8)
