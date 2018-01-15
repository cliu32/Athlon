# Athlon
Accurate Typing of Human Leukocyte Antigen (HLA) by Oxford Nanopore Sequencing
* The pipeline is validated to perform HLA typing using nanopore sequencing reads at 2-field resolution (limiting to the exons 2 and 3 of class I HLA genes at this point). 
* Locus-specific sequencing reads should cover at least exons 2 and 3 of HLA-A, -B, and -C genes, and the number of reads is ideally above 1000 per locus.  

## Prerequisites

The following software or compatible versions must be available in the $PATH:
* seqtk, version 1.0-r31
* BLASR, version 2.0.0
* Samtools, version 1.2 using htslib 1.2.1
* Bedtools, version 2.25.0
* Pyfaidx, version 0.4.8.3
* freebayes, version 1.1.0
* vcflib, version 1.0.0
* Blast, version 2.2.31

## Getting Started and analyzing the test samples

```
git clone https://github.com/cliu32/Athlon.git
cd Athlon
chmod +x athlon.sh
./athlon.sh <ReadNumber>
```

Three samples in the /data folder will be analyzed, each containing ~3000 reads for locus A, B, and C, respectively. For <ReadNumber> when calling athlon.sh, please specify the number of reads (3000 or lower) to be analyzed, for example, 1000. A new folder /data1000 will be generated to hold fastq files with randomly sampled reads. Intermediate files during the analysis will be located the /rslt1000 folder. A .log file and a .rslt file will be generated in the end. The log file will outline the process of selecting candidate alleles; the rslt file will report all the typing results at 3-field G-group level.

## Demultiplexing by locus-specific primer sequences

If you are starting from a fastq file, <fqName>.fastq, that have not been demultiplexed, please follow these steps. i) Include the fastq file in the /rawdata folder, ii) update the primer sequences in the hla_primers_ABC.csv file in the /primer folder to match your primer sequences, and iii) run the following python script: 
```
python kmer_split.py <fqName> hla_primers_ABC
```
Then copy all the demultiplexed fastq files to /data, which will be ready for athlon analysis. 
To customize the demultiplexing step, please follow the help message: 
```
python kmer_split.py -h
```



## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
