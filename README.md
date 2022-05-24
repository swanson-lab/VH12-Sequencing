# VH12-Sequencing
Code and scripts used for Worth et al. 2021: "Receptor editing constrains development of  phosphatidyl choline-specific B cells in VH12-transgenic mice". 

DATA PROCESSING PIPELINE

pRESTO
The data processing pipeline was run in an 8 core Red Hat Enterprise Linux Server with 32Gb RAM was used to process the data.   Upon completing the sequencing run, the entire run folder was transferred from the MiSeq on a portable hard drive to a computer running macOS software from Mojave version 10.14.6 through Catalina version 10.15.3.  The sample sheet for the run was edited and saved as a .txt file in Microsoft Excel (version 16.35) to include the correct indices and sample names for accurate demultiplexing, along with specifying the Read 1 and Read 2 adapters for adapter trimming. FileZilla (versions 3.43.0 through 3.47.2.1 for macOS) was used to perform SFTP (Secure File Transfer Protocol) to upload the entire sequencing data folder to a run-specific directory on the server. Subdirectories were created within the run directory for each sample to separate duplicitously named files generated during the analysis of each sample data set. The code (available at https://github.com/swanson-lab/VH12-Sequencing/) was executed by both copying sections directly into the terminal and being run as a bash script. 
After the raw MiSeq run folder with the updated Sample Sheet was uploaded and entered, then the bcl2fastq program (v2.20, Illumina) was used with parameters `"--use-bases-mask Y*,I6N2,N2I6,Y* --barcode-mismatch 0"` to convert the raw read data into FASTQ files, demultiplex samples based on the indices provided in the sample sheet, and remove the Illumina adapters used to bind the sequence amplicons to the flow cell.
Each set of paired end reads was quality filtered using FilterSeq.py from pRESTO (v0.5.13, Immcantation) with parameters `"-q 20"`. For each sample, the V region primers were masked on read 1 and the C region primers were cut on read 2, both using MaskPrimers.py, with the respective parameters `"-p ../V_primers.fasta --start 0 --mode mask"` and `"-p ../C_primers.fasta --start 15 --mode cut --barcode."` Reads were then subjected to the first round of sequence pairing to sort and match sequence records with equivalent coordinates across files using PairSeq.py with parameters `"--2f BARCODE --coord illumina."` After pairing, a consensus sequence was determined for each set of reads from the sample using BuildConsensus.py with parameters `"--bf BARCODE --pf PRIMER --prcons 0.6 --maxerror 0.1 --maxgap 0.5"` for read 1 and `"--bf BARCODE --pf PRIMER --maxerror 0.1 --maxgap 0.5"` for read 2. Reads were subjected to a second round of sequence pairing using PairSeq.py to again sort and match sequence records with the same coordinates across files with parameters `"--coord presto --1f CONSCOUNT --2f CONSCOUNT PRCONS."` The reads were then aligned and condensed into a single sequence; this required that read 2 be transformed into the reverse complement during the process and was performed using AssemblePairs.py with the parameters `"align --coord presto --rc tail --1f CONSCOUNT --2f CONSCOUNT PRCONS."` Next, the pRESTO annotations were parsed into FASTA/FASTQ sequence headers using ParseHeaders.py with parameters `"collapse -f CONSCOUNT --act min."` Duplicate sequences from FASTA/FASTQ files were removed using CollapseSeq.py with parameters `"-n 20 --inner --uf PRCONS --cf CONSCOUNT --act sum."` The pRESTO annotations were again parsed into FASTA/FASTQ sequence headers using ParseHeaders.py with parameters `"table -f ID PRCONS CONSCOUNT DUPCOUNT."` The records in the console log of the pRESTO modules were then parsed using ParseLog.py for each set of logs with the respective parameters `"-l FS1.log FS2.log -f ID QUALITY"`,  `"-l MP1.log MP2.log -f ID PRIMER BARCODE ERROR"`, `"-l BC1.log BC2.log -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR"`, and `"-l AP.log -f ID LENGTH OVERLAP ERROR PVALUE FIELDS1 FIELDS2."` The FASTQ sequence files were then sorted, sampled, and split into groups with 2+ or one UMI using SplitSeq.py with the parameters `"group -f CONSCOUNT --num 2 --fasta."`

IMGT/HighV-Quest
After pairs of FASTA files for each sample were generated, the FASTA file for sequences with 2 or more copies mapping back to a single UMI (identified with the suffix "_atleast-2" in the file name) was uploaded to IMGT/HighV-Quest using the following run parameters: `"Analysis title - "SampleN 2," Species - Mus musculus (house mouse), Receptor type or locus - IG, Upload sequences in FASTA format - "SampleN_atleast-2.fasta.""` One analysis run was performed for each individual sample.

ChangeO
Once the IMGT/HighV-Quest run was completed (this normally took 6 hours but could take at least 24 hours), the .txz file for each sample was downloaded from the analysis page. These were then uploaded to the linux server in a new folder for the sequencing run along with the original FASTA files submitted for the analysis. Work was again performed in the terminal shell and began with downloading up-to-date germlines from IMGT for reference during processing by entering `"fetch_imgtdb.sh -o ~share/germlines/imgt."` The Change-O package (v0.4.6, Immcantation) was utilized to generate a data table from both the IMGT and FASTA files for each sample using MakeDb.py with parameters `"imgt --regions --scores."` Functional sequences were selected for using ParseDb.py with parameters `"select -f FUNCTIONAL -u T."` Finally, the full germline sequences were added to the data table for comparison using CreateGermlines.py with parameters `"-g full -r ~/share/germlines/imgt/mouse/vdj/*.fasta."` While preparing this manuscript it was noted that ChangeO v1.1.0 defaults to create data tables formatted based on AIRR (Adaptive Immune Receptor Repertoire) guidelines and thus differ slightly in style but not content from the data tables referenced in this work.

Alakazam for V and J gene usage
In order to complete the data analysis and visualization steps, work was moved to RStudioServer running R (v 3.6.0). For each replicated data set, a new R project was created and saved in the folder that the ChangeO work was performed in. The functions gsubGene2NumV, gsubGene2NumJ, gsubNum2GeneV, and gsubNum2GeneJ were copied and saved as objects in the R environment. Libraries were loaded for each package used. The Alakazam package (0.3.0, Immcantation) was utilized to generate an R dataframe from the ChangeO files using function readChangeoDb(). The usage frequencies of the V segment genes were tallied using `countGenes(data, gene = "V_CALL", mode = "gene").` The type of sample was added to the data table for identification of genotype and phenotype using the dplyr (0.8.3) function `mutate(data, Sample = "___").` The list of genes was reordered based on their location within the locus according to IMGT using the gtools (3.8.1) function `mixedorder(data$GENE, decreasing=TRUE)` and the `gsubGene2NumV(data$GENE)` and `gsubNum2GeneV(data$GENE)` functions. Factors and levels were then reassigned to maintain the preferred order by specifying `factor(data$GENE, levels=data$GENE)`. Representation of all genes on the light chain locus was ensured for each sample using stringr (1.4.0) to detect strings of interest and tibble (2.1.3) to modify the data table, along with the V_GENE vector in order to create graphs with axes containing the full Vk and Vl loci. A new data frame for each sample type was created to calculate the mean usage frequency and standard error of the mean (SEM) of each gene in the repertoire across the biological replicates for each sample. First the sequence frequency column for each sample was selected, then the values were averaged and used to determine the SEM before the data frame was annotated with the respective light chain V gene calls and the appropriate sample types. A single, summative data frame with the mean frequency, SEM, gene names, and sample names was then generated by combining the individual sample type data frames. The list of genes was reordered based on their location within the locus according to IMGT and factors and levels for the genes were reassigned as performed above. Finally, factors and levels of the sample types were manually reassigned to maintain the desired specified order within the graph. Four commonly used pseudogenes were labeled with ψ.
Using ggplot2 (3.2.1), a graph of total gene usage by sample type was drawn with the V gene call on the x axis and mean frequency±SEM in the repertoire on the y axis, applying `geom_col()` and `geom_errorbar()` with `position_dodge()` for spacing. The process was repeated for J gene usage and the respective graph was generated.

V gene specific CDR3 and associated J gene usage analysis
Work was continued in RStudioServer with the data tables created from the ChangeO processing for the addition the amino acid sequence for the CDR3 region using `mutate()` from dplyr and `translateDNA()` from Alakazam. Sub-tables for each sample were generated based on the specified V gene call of interest, starting with IGKV4-91, by using `str_detect(data$V_CALL, "4-91").` At this point, only necessary columns were pulled with select(data, column names). The usage frequencies of J gene segments among sequences using IGKV4-91 were tallied using `countGenes(data, gene = "J_CALL", mode = "gene").` As above, `mixedorder()` and the `gsubGene2NumJ()` and `gsubNum2GeneJ()` functions were applied to set the correct order of the genes based on position in the locus before factors and levels were reassigned to ensure maintenance of the new order. Verification that all J segment genes on the kappa locus were represented for each sample, to create graphs with identical x axes, occurred by adding any missing segments with usage frequencies of 0.0 as needed. A new data frame was generated to calculate the mean usage frequency and standard error of the mean for each gene in the J kappa locus across the biological replicates for all sequences using IGKV4-91 from each sample by first selecting the sequence frequency column from individual replicate tables, then averaging the values and calculating the SEM before annotating the data frame with the respective light chain V gene calls and appropriate sample types. A single summative data frame with the mean frequency, SEM, gene names, and sample names was assembled by combining the individual sample type data frames. The list of genes was reordered based on their location within the locus according to IMGT and factors and levels for the genes were reassigned as before. Finally, manual reassignment factors and levels of the sample types was performed to maintain the desired specified order of the samples within the graph. Again using ggplot2, a graph of total gene usage by sample type was drawn with J gene call on the x axis and mean frequency in the repertoire ± SEM on the y axis by applying `geom_col()` and `geom_errorbar()` with `position_dodge()` for spacing.
 
Sequence logos for the most common length of the CDR3 for each subset of sequences using one of the Vk genes of interest (i.e. IGKV4-91) were generated using ggseqlogo (0.1).  The mode CDR3 length for sequences using IGKV4-91 was determined using the function getmode() and sequences with that CDR3 length were selected for using `mutate()`, `str_length()`, and `filter()`. One data frame was created for each sample type with all the CDR3 sequences of the mode length from sequences using IGKV4-91 from the three biological replicates. A sequence logo with the mode CDR3 length was then drawn for each sample type from the sequences using IGKV4-91 with `ggseqlogo(data$CDR3_AA, method = 'prob')`. Ggseqlogo is fully compatible with ggplot2 so additional modifications to the charts were made with the theme functions from ggplot2. This process was repeated for sequences using IGKV3-7, IGKV4-86, IGKV1-110, and IGKV14-126. For the samples comparing the bone marrow repertoire to that of the splenic populations, this analysis was also repeated for sequences using IGKV5-43, and IGKV1-135.

Repertoire Dissimilarity Index (RDI)
All repertoires analyzed in the study were compared by calculating the Repertoire Dissimilarity Indices among pairs of repertoires, using the RDI package in Immcantation (1.0.0). Gene usage information from all samples was combined into one table, with the genotype, phenotype, and biological replicate included as annotations using `mutate()`. The RDI was determined by first defining variables comprising the genes in the repertoire and the sequence annotations describing the replicate and sample type. The gene counts and frequencies were tabulated using `calcVDJcounts()`. This output was then used to generate the RDI matrix through `calcRDI(data, distMethod = "euclidean")`.  A dendrogram showing how much dissimilarity there is among the various repertoires was drawn using plot(hclust(data)). A heat map highlighting the log fold change differences among the repertoires from the RDI values was constructed in Microsoft Excel (v16.37).

 
LIST OF FILES, HARDWARE, AND SOFTWARE

Annotated Code
* Annotated Sample pRESTO Code
* Annotated Sample Change-O Code
* Annotated Assembly of Combined Graphs
* Annotated Assembly of Combine Vk Graphs
* Annotated Combined Sequence Logos
* Annotated NP and CDR3 Length Anovas
* Annotated RDI Code

Sample Files (All data can be found at GEO: https://www.ncbi.nlm.nih.gov/geo/. Accession number: GSE165776)
* Sample Sheet Example
* Sample pRESTO output WT_5N_LN.fasta
* Sample IMGT output WT_5N_LN.txz
* Sample annotated functional sequences WT_5N_LN.tab.zip

Functions and Primers
* C_primers.fasta
* V_primers.fasta
* V_GENES.txt
* J_GENES.txt
* getmode.R
* sem.R
* pseudogenes.R
* gsubGene2NumJ.R
* gsubGene2NumV.R
* gsubNum2GeneJ.R
* gsubNum2GeneV.R

Hardware Used
* Illumina MiSeq™ Instrument ID M03509, CPU ID HWI-M03509 and Control Software version 2.6.2.1
* Agilent 2100 Bioanalyzer with software version B.02.09.SI725 (SR1); Qubit Fluorimeter 3.0 SN 2321602292
* Thermo Scientific Arktik Thermal Cycler V01.07, 48-dual block, SN AKC011203807; 8 core (each on Intel® Xeon® CPU E5-2680 v4 @ 2.40GHz) 
* Red Hat Enterprise Linux Server (release number 3.10.0-514.6.2.el7.x86_64; version #1 SMP Thu Feb 23 03:04:39 UTC 2017, hardware architecture x86_64 [64bit], operating system GNU/Linux) with 32Gb RAM, set up by the Creighton University RadLab.

Software Used
* FileZilla version 3.43.0 for macOS (application)
* Python 2.7.5 (default, Jun 20 2019, 20:27:34) [GCC 4.8.5 20150623 (Red Hat 4.8.5-36)] on linux2 (language)
* R version 3.6.0 (2019-04-26) -- "Planting of a Tree" (language), Copyright (C) 2019 The R Foundation for Statistical Computing
* R Studio Server Swanson Computer port 8787, v1.2.5001 (server application)
  * Stringr version 1.4.0
  * Dyplr version 0.8.3
  * Ggplot2 version 3.2.1
  * Ggseqlogo version 0.1
  * Gtools version 3.8.1
  * Scales version 1.1.0
* bcl2fastq v2.20.0.422 Copyright (c) 2007-2017 Illumina, Inc.
* Immcantation Programs
  * pRESTO version 0.5.13 (program on server)
  * Changeo version 0.4.5 (program on server)
  * Alakazam version 0.3.0 (program on server)
  * RDI version 1.0.0 (program on server)
* IMGT/HighV-Quest version 1.6.7 (web-based IMGT portal)
* Devtools version 2.2.1
* Xml 2 version 1.2.2
* Dependencies for Immcantation programs
  * Python version 3.7.2
  * Biopython version 1.7.3
  * setuptools version 2.0
  * NumPy version 1.8
  * SciPy version 0.14
  * pandas version 0.15
  * airr version 1.2.1
  * Phylip version 3.697




