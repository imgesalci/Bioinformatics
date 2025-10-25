# 🧬 Protein Structure Analysis using R – Insulin (9INS) & Keratin (3ASW)

This repository contains my bioinformatics project analysing the structures of **Insulin** and **Keratin** proteins using **R**, and **Bio3D**.  
The aim was to understand their amino acid composition, identify similar structures with **BLAST**, perform **Normal Mode Analysis (NMA)**, and visualize them in 3D.

---

## 🔍 Project Overview

This project explores:
- **Amino acid composition** using sequence data from UniProt.  
- **Structural similarity** via BLAST searches on RCSB PDB.  
- **Normal Mode Analysis (NMA)** to investigate flexibility and molecular motion.  
- **3D visualization** using the `rgl` package to display spatial organization of atoms.

All analysis was performed using **R**.

---

## 🧩 Methodology

### 1️⃣ Data Retrieval
Protein sequences and structures were retrieved from:
- [UniProt – P01308 (Insulin)](https://www.uniprot.org/uniprot/P01308)
- [UniProt – P04264 (Keratin)](https://www.uniprot.org/uniprot/P04264)
- [RCSB PDB – 9INS (Insulin)](https://www.rcsb.org/structure/9INS)
- [RCSB PDB – 3ASW (Keratin)](https://www.rcsb.org/structure/3ASW)

```r
library(bio3d)
library(seqinr)

insulin.seq <- read.fasta("https://www.uniprot.org/uniprot/P01308.fasta")
insulin.pdb <- read.pdb("https://files.rcsb.org/view/9INS.pdb")

keratin.seq <- read.fasta("https://www.uniprot.org/uniprot/P04264.fasta")
keratin.pdb <- read.pdb("https://files.rcsb.org/view/3ASW.pdb")
```

---

### 2️⃣ Amino Acid Composition
Calculated the amino acid frequencies for both proteins and visualized them using bar plots.
 
<img width="640" height="480" alt="compositions" src="https://github.com/user-attachments/assets/7ba993b9-5833-4c8e-9c7f-6c8fd7b26d76" />


---

### 3️⃣ BLAST Similarity Search
Performed BLAST to find the 10 most similar structures for each protein using:

```r
insulin.blast <- blast.pdb(insulin.pdb)
keratin.blast <- blast.pdb(keratin.pdb)
plot.blast(insulin.blast, top = 10)<img width="1450" height="928" alt="insulin1" src="https://github.com/user-attachments/assets/75b612c3-938d-4795-ba58-5b43b037101c" />

plot.blast(keratin.blast, top = 10)
```

📊 *Blast Results:*  
<p align="center"> <img width="480" height="360" alt="insulin blast results" src="https://github.com/user-attachments/assets/397635e0-4722-4a5b-bf83-6b4c032ef7dd" /> <br> <em>Figure 3. BLAST similarity results for Insulin (9INS)</em> </p> <p align="center"> <img width="480" height="360" alt="keratin blast results" src="https://github.com/user-attachments/assets/237488f4-2f05-4254-9ec8-bfa27113583e" /> <br> <em>Figure 4. BLAST similarity results for Keratin (3ASW)</em> </p>
---

### 4️⃣ Normal Mode Analysis (NMA)
Analyzed molecular motion and flexibility.

```r
insulin.modes <- nma(insulin_trimmed)
keratin.modes <- nma(keratin_trimmed)
plot.nma(insulin.modes)
plot.nma(keratin.modes)
```

🧠 *Insight:* Insulin showed higher eigenvalues (more rigid structure), while keratin displayed greater flexibility, especially at molecular ends.

---

### 5️⃣ 3D Visualization
Visualized both structures using `rgl` for spatial representation.

```r
open3d()
plot3d(insulin_coords$x, insulin_coords$y, insulin_coords$z, col = "blue", size = 1)
title3d("3D Visualization of Insulin")
```

🧩 *3D Images:*  
<p align="center"> <img width="480" height="360" alt="insulin 3D visualization" src="https://github.com/user-attachments/assets/0c2fc4cf-4608-4baf-82cb-744e6f138213" /> <br> <em>Figure 5. 3D spatial model of Insulin</em> </p> <p align="center"> <img width="480" height="360" alt="keratin 3D visualization" src="https://github.com/user-attachments/assets/c183afc9-6345-47a7-acc9-a8007d2c0120" /> <br> <em>Figure 6. 3D spatial model of Keratin</em> </p>


---

## ⚙️ Tools & Libraries
- **R packages:** bio3d, seqinr, rgl
- **Databases:** UniProt, RCSB PDB  
- **Techniques:** BLAST, NMA, 3D rendering  

---

## 🎯 Key Findings
- Insulin and keratin have distinct amino acid distributions.  
- Insulin’s structure is more compact and rigid.  
- Keratin exhibits higher flexibility at terminal regions.  
- Both share sequence similarities across species, suggesting evolutionary relationships.

---

📌 *Developed by **Elif İmge Şalcı** | Ankara University, 2024*
