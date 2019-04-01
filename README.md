# LG_to_Chromosome_SyntenyPlots
This repository contains a custom R script to plot Circos plots using Circlize to identify synteny between linkage groups and chromosomes in a reference genome

three input files are needed: 
 1. FalconChromosomeLenghts.txt - a list of chromosome lengths 
       - in my case this included the lengths of the fasta headings (all < 50) followed by the number of bases in the chromosome fasta line

e.g. 
45
111295376
45
110734514
45
108247109
44
92930148

 2. LM_stats.txt - a summary statstics table of the linkage (taken from Table 1: http://www.g3journal.org/content/8/12/3745.figures-only)

e.g. 
Calb01 253 75.96 0.30 91.07 0.36 63.67 0.25 Ssa01 W02 1.43
Calb02 228 83.57 0.37 101.33 0.44 69.58 0.31 Ssa01 W03 1.46
Calb03 220 78.51 0.36 84.40 0.38 87.95 0.40 Ssa21 W32 0.96

 3. SA_linkagemap_FalconPhaseMapped_Filtered.csv - this is the result of mapping all linkage map rad loci to the reference and filtering for MAPQ > 30 using: 
bsub -n10 -W 4:00 -R "rusage[mem=3000]" "module load java gdc bwa/0.7.17; bwa mem -t 10 module load java gdc bwa/0.7.17; bwa mem -t 10 chromosome_level_assembly.fasta SpeciesAveragedLinkageMap.fasta > SA_linkagemap_FalconPhaseMapped.sam"

awk '{ if($5 > 30) {print}}' SA_linkagemap_FalconPhaseMapped.sam > SA_linkagemap_FalconPhaseMapped_Filtered.txt

this was then checked over in excel and the column with mapping location separated by scaffold number and bp location 
       e.g. check columns carefully first, there are some gaps in the header line
@PG			ID:bwa	PN:bwa		VN:0.7.17-r1188	CL:bwa	mem#NAME?	10	/cluster/project/gdc/shared/p298/Rishi/PostPhaseScaffoldingGenomes/WF_Canu.chr.fasta	/cluster/project/gdc/shared/p298/Rishi/LinkageMap/SpeciesAveragedMap.fasta						
consensus_42108	SA10	0	0	PGA	scaffold21	16735146	60	90M	*	0	0	CACTCTCTCTCCTTCTCCCTCCTTTCCTCTCTATGTCCCAGAAAGCCTTCAAAGATGACCCGCAGAGCAGAGCGGTGGTATGATCCTGCA	*	NM:i:0	MD:Z:90	AS:i:90	XS:i:21	
consensus_27024	SA10	2.3	0	PGA	scaffold21	6158036	60	90M	*	0	0	TGCAGGGTTCAAATGCAACCACTTGGGCAGCCTGGTCTCATAGGCTTGACGTAACATAGAAAATGTAAATCTGGGACACTCAAATTAGTA	*	NM:i:0	MD:Z:90	AS:i:90	XS:i:23	
consensus_7502	SA10	2.63	16	PGA	scaffold21	26075329	60	90M	*	0	0	GCAGCTGCTGCTTGTGGGCCTCCTTAATTGCTTGCAGCTGGCGTTGGAGCAGCAGCTCTGATCTCAGGCCATCCTCCAGGAAAGCCTGCA	*	NM:i:1	MD:Z:25C64	AS:i:85	XS:i:20	

for extra information about plotting check out the circlize manual: https://jokergoo.github.io/circlize_book/book/ 

