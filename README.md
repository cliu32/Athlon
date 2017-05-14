# Athlon
Accurate Typing of Human Leukocyte Antigen (HLA) by Oxford Nanopore Sequencing
* The pipeline is validated to perform HLA typing using nanopore sequencing reads at 3-field resolution with the G-group nomenclature. 
* Locus-specific sequencing reads should cover at least exons 2 and 3 of HLA-A, -B, and -C genes, and the number of reads is ideally above 100 per sample. 
* Better results have been achieved with R9.4 flow cell chemistry than R7.4 flow cell chemistry. 

## Prerequisites

The following software or compatible versions must be available in the $PATH:
* BLASR, version 2.0.0
* Samtools (Samtools and bcftools), version 1.2 using htslib 1.2.1
* Bedtools, version 2.25.0
* Pyfaidx, version 0.4.8.3
* Blast, version 2.2.31

## Getting Started and running the tests

```
cd Athlon
chmod +x athlon.sh
./athlon.sh
```

Six samples in the data folder will be analyzed. A .log file and a .rslt file will be generated. Intermediate files will be available in the rslt folder upon completion of the analysis. 


## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
