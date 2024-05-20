# envs
The tools needed to used in shell scripts.

# Shell scripts
Shell scripts used in genome assembly, viral recognition and other analysis.

# Rawdata processing
Raw next generation sequencing of viral reads were processed with **Trimmomatic v0.38** (with parameter LEADING:3 TRAILING:3 SLIDINGWINDOW:15:30 MINLEN:50) to remove adaptors and trim low-quality bases; reads of 50 bp or less after trimming were discarded. Putative human reads were identified from the trimmed reads by aligning the latter to the human reference genome (hg38; GCA_000001405.15) using **Bowtie2** (v2.4.2, --end-to-end) with default parameters and removed from further analysis.

# Assembly
* There is two method used in our study, The main method in manuscript is consistent with Zeng S, et al. The quality-controlled reads of bulk and VLP were assembled with **MegaHIT v1.1.3** (default parameters except option “-min-contig-len 1000”).
#
* For our method (Hereafter is our_method), Before assembly, Subsequently, we removed the bacterial contamination by aligning both bulk and VLP clean reads to UHGG-Minus genomes (the Unified Human Gastrointestinal Genome (UHGG) catalog that removed possible prophage regions) using **Bowtie2 v2.4.2**. **IDBA-UD** (**Release 1.1.3**, parameters: --maxk 120 --step 10 –min_contig 1000) was used to assemble the filtered bulk and VLP data in each sample. 
# Viral identification and acquisition of non-redundant vOTUs
* The method for viral gene identification is consistent with Zeng S, et al. After assembly, all contigs > 3kb in length in each sample were submitted to **VirSorter v2.2.3, VirFinder v1.1 and VIBRANT v1.2.1** for identification of viral populations (The parameters is the same as the Methods in Zeng S, et al.7). Thereafter, the viral populations were filtered by **CheckV v1.0.1**, and further dereplication by **CD-HIT v4.6.8** (parameters: -c 0.95 -M 16000 -aS 0.85 -d 0) to obtain non-redundant species-level viral operational taxonomic units (vOTUs) based on the two methods. 
#
* For our_method, we use another method for viral gene identification is referenced in the GPD database. After assembly, all contigs > 5kb in length in each sample were submitted to **VirSorter v2.0** (--min-score 0.7) and **VirFinder v1.1** (default parameters) for identification of potential viral operational taxonomic units (vOTUs). We consider a contig as a virus if it satisfies both classification criteria. Next, we merged all the vOTUs from all samples in mNGS and vNGS separately, and then the merged dataset was dereplicated using **CD-HIT** (**v4.8.1**, parameters: -c 0.95 -M 16000 -aS 0.8 -d 0) to establish two non-redundant viral genome catalogs based on the two methods. 

# R scripts
Shell scripts used in abundance calculation and figure generation.
