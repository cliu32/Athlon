#Athlon: Accurate typing of human leukocyte antigen (HLA) by Oxford Nanopore sequencing
#Author: Chang Liu, Washington University School of Medicine, cliu32@wuslt.edu
#!/bin/bash
{
start=$(date)
datafolder=data/
rsltfolder=rslt/ #will hold intermediate files after completion of the program.
reffolder=reference/

mkdir $rsltfolder
> $rsltfolder'rslt'.txt #will be moved to the pwd after completion of the program.

for sample in $(ls $datafolder | awk -F"[_]" '{print $1"_"$2}')
do
    seqfile=$datafolder$sample'*.fastq'
    locus=$(echo $sample | cut -d'_' -f 2)
    ref='g_3f_hla3260_'$locus'_exon_2_3_introngap'
    
    echo '*********||| sample: '$sample', locus: '$locus' |||*********' 
    echo $sample >> $rsltfolder'rslt'.txt 
    echo 'indexing reference '$reffolder$ref'.fasta ...' 
    samtools faidx $reffolder$ref'.fasta' 
    
    echo 'aligning reads to the reference...' 
    blasr $seqfile $reffolder$ref'.fasta' -hitPolicy randombest -sam -clipping soft -out $rsltfolder$sample.sam 
    
    echo 'converting sam file to bam file, sorting, and indexing ...' 
    samtools view -b $rsltfolder$sample.sam -o $rsltfolder$sample.bam
    samtools sort $rsltfolder$sample.bam $rsltfolder$sample.sorted
    samtools index -b $rsltfolder$sample.sorted.bam
    
    echo 'quantify coverage by bedtools ...' 
    echo 'generating bed file from the reference file ...' 
    faidx --transform bed $reffolder$ref.fasta > $reffolder$ref.bed 
    echo 'sum coverage of 3-field typing and sort in descending order ...' 
    coverageBed -abam $rsltfolder$sample.sorted.bam -b $reffolder$ref.bed -d | groupBy -g 1 -c 14 -o sum | sort -rV -k 2,2 > $rsltfolder$sample$'_3f_coverage'.csv
    
    echo 'sum coverage of 2-field typing and sort in descending order ...' 
    awk -F":" '$1=$1":"$2' OFS=" \t" $rsltfolder$sample$'_3f_coverage'.csv | sort -k 1,1 | groupBy -g 1 -c 4 -o sum | sort -Vr -k 2,2 > $rsltfolder$sample$'_2f_coverage'.csv
    
    echo 'sum coverage of 1-field typing and sort in descending order ...' 
    awk -F":" '$1=$1' OFS="\t" $rsltfolder$sample$'_3f_coverage'.csv | sort -k 1,1 | groupBy -g 1 -c 4 -o sum | sort -Vr -k 2,2 > $rsltfolder$sample$'_1f_coverage'.csv
    
    echo '***top 10 1-field typings***' 
    head $rsltfolder$sample$'_1f_coverage'.csv 
    
    echo '***top 10 2-field typings***' 
    head $rsltfolder$sample$'_2f_coverage'.csv 
    
    echo '***top 10 3-field typings***' 
    head $rsltfolder$sample$'_3f_coverage'.csv 

    echo 'generating 1-field candidate allele pair ...' 
    allele0=$(awk 'NR==1 {print $1}' $rsltfolder$sample$'_1f_coverage'.csv) #1-field dominant allele
    cover0=$(awk 'NR==1 {print $2}' $rsltfolder$sample$'_1f_coverage'.csv)
    allele1=$(awk 'NR==2 {print $1}' $rsltfolder$sample$'_1f_coverage'.csv) #1-field minor allele
    cover1=$(awk 'NR==2 {print $2}' $rsltfolder$sample$'_1f_coverage'.csv)
    
    echo $allele0 
    if (( $(( 100*cover1 )) < $(( 15*cover0 )) )) #if the 1-field minor allele has < 15% coverage of the 1-field dominant allele, call the 1-field dominant allele only
    then 
        echo '-'  
        echo 'generating 2-field candidate allele pair ...' 
        cat $rsltfolder$sample$'_2f_coverage'.csv | awk -v a="$allele0" -F":" '$1==a' | head -n 2 > $rsltfolder$sample'candi_2f'.csv
        allele0=$(awk 'NR==1 {print $1}' $rsltfolder$sample'candi_2f'.csv) #2-field dominant allele
        cover0=$(awk 'NR==1 {print $2}' $rsltfolder$sample'candi_2f'.csv)
        allele1=$(awk 'NR==2 {print $1}' $rsltfolder$sample'candi_2f'.csv) #2-field minor allele
        cover1=$(awk 'NR==2 {print $2}' $rsltfolder$sample'candi_2f'.csv)
        
        echo $allele0 
        if (( $(( 10*cover1 )) < $(( 4*cover0 )) )) #if the 2-field minor allele has < 40% coverage of the 2-field dominant allele, call the 2-field dominant allele only
        then 
            echo '-'  #'there seems to be only one 2-field candidate allele: '
            sed -i /$allele1/d $rsltfolder$sample'candi_2f'.csv #remove the 2-field minor allele 
        else
            #do nothing, keep the 2-field minor allele
            echo $allele1 
        fi
    else
        echo $allele1  #print the 1-field minor allele
        echo 'generating 2-field candidate allele pair ...' 
        cat $rsltfolder$sample$'_2f_coverage'.csv | awk -v a=$allele0 -F":" '$1==a' | head -n 1 > $rsltfolder$sample'candi_2f'.csv
        cat $rsltfolder$sample$'_2f_coverage'.csv | awk -v a=$allele1 -F":" '$1==a' | head -n 1 >> $rsltfolder$sample'candi_2f'.csv
        cut -f1 $rsltfolder$sample'candi_2f'.csv 
    fi
  
    echo 'generating 3-field candidate allele pair ...' 
    > $rsltfolder$sample'candi_3f.csv' #create an empty file to store 3f candidates
    for allele in `cut -f1 $rsltfolder$sample'candi_2f'.csv`
    do
        grep $(echo "$allele" | sed "s/*/\\\*/g") -m 1 $rsltfolder$sample$'_3f_coverage'.csv | cut -f1 >> $rsltfolder$sample'candi_3f'.csv #Error prone, e.g. searching for A*02:01 results in A*24:02:01G; this is solved by the sed command to replace * with \*
    done
    cat $rsltfolder$sample'candi_3f'.csv 
    
    echo 'preparing fasta files for candidate allele pair ...' 
    > $rsltfolder$sample'candi_3f'.fasta
    for allele in `cat $rsltfolder$sample'candi_3f'.csv`
    do
        samtools faidx $reffolder'g_3f_hla3260_'$locus.fasta $allele >> $rsltfolder$sample'candi_3f'.fasta
    done
    samtools faidx $rsltfolder$sample'candi_3f'.fasta
    
    echo 'realign reads to 3-field candidate allele pair ...' 
    blasr $seqfile $rsltfolder$sample'candi_3f'.fasta -hitPolicy randombest -sam -clipping soft -out $rsltfolder$sample'_3f'.sam
    echo 'converting sam file to bam file, sorting, and indexing ...' 
    samtools view -b $rsltfolder$sample'_3f'.sam -o $rsltfolder$sample'_3f'.bam
    samtools sort $rsltfolder$sample'_3f'.bam $rsltfolder$sample'_3f'.sorted
    samtools index -b $rsltfolder$sample'_3f'.sorted.bam
    
    echo 'generating consensus sequence as a multi-fasta file ...' 
    samtools mpileup -uf $rsltfolder$sample'candi_3f'.fasta $rsltfolder$sample'_3f'.sorted.bam -d 8000 | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 8000 > $rsltfolder$sample'_cns'.fastq
    seqtk seq -a $rsltfolder$sample'_cns'.fastq > $rsltfolder$sample'_cns'.fasta
    
    echo 'blastn the consensus sequence within the reference allelesca ...' 
    cd $rsltfolder
    rm allele*.fasta #clean up the consensus fasta files from the previous sample
    cat $sample'_cns'.fasta | awk '{if (substr($0,1,1)==">") {filename=("allele"x++".fasta")} print $0 > filename}'
    cd ..
    for allele in `ls $rsltfolder'allele'*.fasta`
    do 
        blastn -query $allele -subject $reffolder'g_3f_hla3260_'$locus.fasta -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge | cat | awk '{print $2}' >> $rsltfolder'rslt'.txt
    done    
    
done
end=$(date)
echo 'program started at '$start 
echo 'program ended at '$end 
} > log
mv log $(date -d "today" +"%Y%m%d%H%M%S").log
mv $rsltfolder'rslt'.txt $(date -d "today" +"%Y%m%d%H%M%S").rslt
