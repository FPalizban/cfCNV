# cfCNV
Our main idea: 
Benchmarking the CNV calling pipeline in the area of cfDNA genomic data
Here, we first start with germline variants and will continue with somatic variants in the next step. and about the data, it is genomic NGS data including both short and long reads. short read sequencing data is the one that obtained from illumina sequencing and the long read sequencing is related to PacBio or Oxford Nanopore.

Real dataset: 
for cfDNA illumina sequencing data: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378

for cfDNA PacBio sequencing data: https://ega-archive.org/studies/EGAS00001006609

for cfDNA Oxford Nanopore sequencing data: https://ega-archive.org/studies/EGAS00001004975

Groundtruth  dataset:
 the groundtruth data like NA12878 is not available for cfDNA data. But there is database of related dataset which can be used for more analysis :
https://generegulation.org/cfdna/


Validation data:
the simulation can be used to validate the benchmarking result. a CNV simulator to manipulate a real cfDNA bam file will be implemented and used for further confirmation
and also there are some read visualizing tools that can help to validate the result of CNV callers. these tools are:
Manual assessment of individual candidate SVs remains critical, with many groups using Integrated Genome Viewer (IGV), or more specialized tools including svviz, SplitThreader, Ribbon , or SVPV
and also using public databases:
Very recently, the most comprehensive resource of small and large CNVs derived from WGS is from the Genome Aggregation Database (gnomAD), which was established using numerous CNV and SV callers.
https://www.internationalgenome.org/human-genome-structural-variation-consortium/


there are few benchmarking articles that can be considered for further investigation:
https://www.nature.com/articles/s41587-020-0538-8
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02636-8


Genomic Tools or Methods:
three different main pipelines will be considered for the above-mentioned datasets


Run Benchmarking Experiments:

the Nextflow will be used to implement the pipelines all together.

Evaluation Metrics:

#################################################################################
tools installation

ichorCNA:
$ git clone https://github.com/broadinstitute/ichorCNA.git
$ R CMD INSTALL ichorCNA 
$ git clone https://github.com/shahcompbio/hmmcopy_utils.git
$ cd hmmcopy_utils
$ cmake . 
$ make
$ cd bin
$ ./readCounter --window 1000000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" ~/Documents/tutoring/IH03_filtered.bam > IH03.wig
$ cd --
$ cd ichorCNA
$ Rscript scripts/runIchorCNA.R --id IH03 --WIG IH03.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 --gcWig inst/extdata/gc_hg19_1000kb.wig --mapWig inst/extdata/map_hg19_1000kb.wig --centromere inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --normalPanel inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds --includeHOMD False --chrs "1" --estimateNormal True --estimatePloidy True --estimateScPrevalence True --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ./

 Spectre 

$ git clone https://github.com/fritzsedlazeck/Spectre.git
$ wget https://github.com/brentp/mosdepth/releases/download/v0.3.4/mosdepth && chmod +x ./mosdepth && ./mosdepth -h
$ ./mosdepth -n --fast-mode --by 1000 IH03_1000 IH03.bam  
$ mv  IH03_1000.* mosdepth-0.3.4/IH03_1000
$ python3 spectre.py  CNVCaller  --bin-size 1000  --coverage mosdepth-0.3.4/IH03_1000/ --sample-id IH03_1000 --output-dir IH03_1000-res/ --reference hs37d5.fa.gz
WisecondorX:

first run : $ 
```
pip install -U git+https://github.com/CenterForMedicalGeneticsGhent/WisecondorX
```

b then install the R packages in R 

the commands to run it:
https://github.com/CenterForMedicalGeneticsGhent/WisecondorX?tab=readme-ov-file



QDNAseq:
It is a r package and can be installed in R:
```
BiocManager::install("QDNAseq")
```

```
BiocManager::install("QDNAseq.hg19")
```




BIC-seq2:

 first download the script for normalization and unzip it;
$ wget https://compbio.med.harvard.edu/BIC-seq/NBICseq-norm_v0.2.4.tar.gz
$  tar xvzf NBICseq-norm_v0.2.4.tar.gz
$ cd NBICseq-norm_v0.2.4/
$ make clean
$ make
then the same for the segmentation part which is for the CNV calling;
$ wget https://compbio.med.harvard.edu/BIC-seq/NBICseq-seg_v0.7.2.tar.gz
$ tar xvzf NBICseq-seg_v0.7.2.tar.gz 
$ cd NBICseq-seg_v0.7.2/
$ make clean
$ make



#############################################################
Investigating the data
 first the data was downloaded from different studies and the bam files were checked for their genome reference, which is hg19 for all.
$samtools view -H M12741R.cutadapter4.blasr.bam |grep '^@SQ'c
and then checking the chromosome length with https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&chromInfoPage=

# View basic statistics of BAM file

$ samtools stats your_file.bam
