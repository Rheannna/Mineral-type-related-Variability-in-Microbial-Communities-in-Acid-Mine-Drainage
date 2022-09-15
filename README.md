# Mineral-type-related Microbiomes and Microbial-driven Biogeochemical Cycling in Acid Mine Drainage

#!/bin/bash

#coding:utf-8

#0 set workpath

        mypath=..../data/AMD_meta

#### move to the directory

        for i in `cat result/design.txt  | tail -n +2 | cut -f 1`;do mv data/${i}/1.Cleandata/*gz data/;done
        
####rename

####QC

#1 

####1.2 FastQC

        mkdir -p temp/fastqc_one
        for i in ${mypath}/data/*.fq
        do
            fastqc ${i}
            mv ${mypath}/data/*fastqc.html ${mypath}/data/*fastqc.zip ${mypath}/temp/fastqc_one
        done

#1.3 multiqc

        multiqc -d temp/fastqc_one/ -o temp/multiqc_one

#1.4 remove adapter

        mkdir -p temp/data_before
        for i in `tail -n +2 result/design.txt | cut -f 1`
        do
            cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 50 -O 5 -j 20 -q 20,20 -o temp/data_before/${i}_1.fastq -p temp/data_before/${i}_2.fastq -f fastq data/${i}.1.fq data/${i}.2.fq
        done

#1.5 remove duplications

        mkdir -p temp/data_cut
        mkdir data
        for i in `tail -n +2 result/design.txt | cut -f 1`
        do
            cd-hit-dup -i temp/data_before/${i}_1.fastq -i2 temp/data_before/${i}_2.fastq -o temp/data_cut/${i}_1.fastq -o2 temp/data_cut/${i}_2.fastq -u 50 -e 2
        done
        for i in ${mypath}/temp/data_cut/*.fastq
        do
            fastqc ${i}
            mv ${mypath}/temp/data_cut/*fastqc.html ${mypath}/temp/data_cut/*fastqc.zip ${mypath}/temp/fastqc_two
        done

        multiqc -d temp/fastqc_two/ -o temp/multiqc_two

#1.6 do again until QC is ok


#seqtk 
        trimfq -L 190 1W1.1.fq > 1W1.1.fastq

#align

#2 SPADES

#2.1  spades

        mkdir -p result/spades
        for i in `tail -n +2 result/design.txt| cut -f 1 
        do
            nohup spades.py --meta -k 21,33,55,77,99,127 -1 temp/data_cut/${i}_1.fastq -2 temp/data_cut/${i}_2.fastq -o result/spades/spades_${i}
        done

#quast 

        conda activate metawrap
        mkdir -p result/quast
        mypath=..../meta_all
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
            quast.py -t 20 ${mypath}/result/spades/spades_${i}/scaffolds.fasta --min-contig 500 -o ${mypath}/result/quast/quast_${i}
    stats.sh ${mypath}/result/spades/spades_${i}/scaffolds.fasta > result/quast/stats_${i} 
        done

#remove sequences < 2000bp

        for i in `tail -n +2 result/design.txt| cut -f 1`
        do  
            cutadapt -m 2000 -o ${mypath}/result/spades/spades_${i}/scaffolds_filter_2000.fasta ${mypath}/result/spades/spades_${i}/scaffolds.fasta
        done

#binning

#3 Metabats2 

#3.0 BBmap


        mkdir -p result/BBmap
        mypath=/share/data1/liushuangjiang/huangye/meta_all
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
            bbwrap.sh in1=${mypath}/data_cut/${i}_1.fastq,singletons.fq in2=${mypath}/data_cut/${i}_2.fastq,null\
            out=${mypath}/result/BBmap/${i}.mapped.sam ref=${mypath}/result/spades/spades_${i}/scaffolds_filter_2000.fasta -Xmx100g qin=33 minid=0.97\ local=t covstats=${mypath}/result/BBmap/${i}.constats.txt covhist=${mypath}/result/BBmap/${i}.covhist.txt \
            basecov=${mypath}/result/BBmap/${i}.basecov.txt bincov=${mypath}/result/BBmap/${i}.bincov.txt

#3.1 sam2bam

#Samtools convert sam files to bam files

        mypath=..../meta_all
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
            samtools view -@ 16 -b -S ${mypath}/result/BBmap/${i}.mapped.sam -o ${mypath}/result/BBmap/${i}.mapped.bam
        done

#3.2 bam2sorted.bam

#samtools

        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
            samtools sort -@ 16 -l 9 -O BAM ${mypath}/result/BBmap/${i}.mapped.bam -o ${mypath}/result/BBmap/${i}.sorted.mapped.bam
        done

#samtools resort bam files to get index

        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
            samtools index -b -@ 16 ${mypath}/result/BBmap/${i}.sorted.mapped.bam ${mypath}/result/BBmap/${i}.sorted.mapped.bam.bai
        done
        samtools index -b -@16 1S01Z10.sorted_mapped.bam 1S01Z10.sorted.mapped.bam.bai
        samtools index -b -@16 1W1.sorted_mapped.bam 1W1.sorted.mapped.bam.bai

#3.3 contig coverage

#input = sorted.bam，use Metabat2 jgi_summarize_bam_contig_depths

        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
            jgi_summarize_bam_contig_depths --outputDepth ${mypath}/result/BBmap/${i}.contig.depth.txt ${mypath}/result/BBmap/${i}.sorted.mapped.bam
        done

#3.4 binning

        mkdir -p result/metabat2
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
                nohup metabat2 -m 2000 --unbinned -t 16 -i ${mypath}/result/spades/spades_${i}/scaffolds_filter_2000.fasta -a ${mypath}/result/BBmap/${i}.contig.depth.txt  -o ${mypath}/result/metabat2/${i}/${i}  -v
        done

#4 bins 

#4.1 checkm

        mkdir -p result/checkm

        for i in `tail -n +2 result/design.txt | cut -f 1`
        do 
            nohup checkm lineage_wf -f result/checkm/${i}_checkm.txt -t 16 -x fa ${mypath}/result/metabat2/${i} ${mypath}/result/checkm/${i}
        done
        for i in `tail -n +2 result/design.txt | cut -f 1`
        do 
            nohup checkm qa ${mypath}/result/checkm/${i}/lineage.ms ${mypath}/result/checkm/${i}
        done

#4.2 Bins selection Completeness>50%，Contamination<5%

        for i in `tail -n +2 result/design.txt | cut -f 1`
        do
        tail -n +4 ${mypath}/result/checkm/${i}_checkm.txt | awk -F '[ ;]+' '{if ($14>=70.00&&$15<=5.00) print $2}' > ${mypath}/result/checkm/${i}_filter.txt
        done
        mkdir -p ${mypath}/result/filter_bin
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
            mkdir ${mypath}/result/filter_bin/${i}
            for BIN in `cat ${mypath}/result/checkm/${i}_filter.txt | cut -f 1`
            do 
                cp ${mypath}/result/metabat2/${i}/${BIN}.fa ${mypath}/result/filter_bin/${i}/${BIN}.fa
            done
        done

#4.2.2 

        for i in `tail -n +2 result/design.txt | cut -f 1`
        do
                tail -n +4 ${mypath}/result/checkm/${i}_checkm.txt | awk -F '[ ;]+' '{if ($14>=70.00&&$15<=5.00) print $2,"\t", $14,"\t",$15}' > ${mypath}/result/checkm/${i}_filter1.txt
        done
        mv result/checkm/*_filter1.txt result/bin_tax/

#4.3 bins coverage

        mkdir -p ${mypath}/result/bin_quant
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
                nohup metawrap quant_bins -t 32  -b ${mypath}/result/filter_bin/${i}  -o ${mypath}/result/bin_quant/${i}  -a ${mypath}/result/spades/spades_${i}/scaffolds_filter_2000.fasta ${mypath}/data_cut/${i}_1.fastq ${mypath}/data_cut/${i}_2.fastq
        done

#####kallisto

        conda activate kallisto

        for i in `tail -n +2 result/design.txt| cut -f 1`
        do
            mkdir -p ${mypath}/result/bin_quant_kallisto/${i}
            for BIN in `cat ${mypath}/result/checkm/${i}_filter.txt | cut -f 1`
            do 
            nohup kallisto index --index=${mypath}/result/bin_quant_kallisto/index/${BIN}.index ${mypath}/result/filter_bin/${i}/${BIN}.fa
            done
        done

        for BIN in `ls ${mypath}/result/bin_quant_kallisto/index/*.index | cut -d "/" -f 11 | cut -d "." -f 1,2`
        do
            mkdir -p ${mypath}/result/bin_quant_kallisto/${BIN}
            for i in `tail -n +2 result/design.txt| cut -f 1`
            do 
            nohup kallisto quant --index=${mypath}/result/bin_quant_kallisto/index/${BIN}.index --output-dir=${mypath}/result/bin_quant_kallisto/${BIN}/${i} --threads=16 --plaintext ${mypath}/data_cut/${i}_1.fastq ${mypath}/data_cut/${i}_2.fastq
            done
        done

#rename

        for BIN in `ls ${mypath}/result/bin_quant_kallisto/index/*.index | cut -d "/" -f 11 | cut -d "." -f 1,2`
        do
            for i in `tail -n +2 result/design.txt| cut -f 1`
            do
            mv ${mypath}/result/bin_quant_kallisto/${BIN}/${i}/abundance.tsv ${mypath}/result/bin_quant_kallisto/${BIN}/${i}.tsv
            done
        done

#combine all bins relative abundance to create mOTU table

        for i in `tail -n +2 result/design.txt| cut -f 1 `;do paste ${mypath}/result/bin_quant/${i}/bin_abundance_table.tab > ${mypath}/result/bin_quant/${i}_bin_abundance_table.tab;done

        cat *bin_abundance_table.tab > all_bin_abundance_table.txt

###################################################    4.4 taxonomy

#gtdbtk1.3

        conda deactivate 
        conda activate gtdbtk1.3
        mkdir -p ${mypath}/result/gtdbtk
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
                gtdbtk classify_wf --cpus 24 --genome_dir ${mypath}/result/filter_bin/${i} --out_dir ${mypath}/result/gtdbtk/${i} --extension fa
        done
        #
        for i in `tail -n +2 result/design.txt| cut -f 1 `;do cat ${mypath}/result/gtdbtk/${i}/*.summary.tsv > ${mypath}/result/bin_tax/${i}.gtdbtk.summary.txt;done
        #
        for i in `tail -n +2 result/design.txt| cut -f 1 `;do awk '!/user_genome/' ${mypath}/result/bin_tax/${i}.gtdbtk.summary.txt > ${mypath}/result/bin_tax/${i}.gtdbtk.summary1.txt;done
        #
        for i in `tail -n +2 result/design.txt| cut -f 1 `;do awk '{ if (NR==FNR) {arraya[$1]=$2} if (NR!=FNR) { arrayb[$1]=$0}}END{for (i in arraya) {print i,'\t',arraya[i],'\t',arrayb[i]}} ' ${mypath}/result/bin_tax/${i}.gtdbtk.summary1.txt ${mypath}/result/bin_tax/${i}_filter1.txt > ${mypath}/result/bin_tax/${i}_gtdbtk_checkm.txt;done

##############################################           4.6 bins rebin

        mkdir -p ${mypath}/result/bin_reassembly
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
          nohup metawrap reassemble_bins -o ${mypath}/result/bin_reassembly/${i} -1 temp/data_cut/${i}_1.fastq -2 temp/data_cut/${i}_2.fastq -t 8 -m 800 -c 50 -x 10 -b ${mypath}/result/filter_bin/${i}
        done
        echo "reassembly done!"

#############################################        4.6.1 remove contamination 
sequences

        contig=${mypath}/result/spades/spades_${i}/scaffolds_filter_2000.fasta
        bins=${mypath}/result/filter_bin/${i}
        bam=${mypath}/result/BBmap/${i}.sorted.mapped.bam
        mkdir -p ${mypath}/result/bin_refinem/scaffold_stats
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do
            mkdir -p ${mypath}/result/bin_refinem/${i}
            #refinem scaffold_stats -c 16 <scaffold_file> <bin_dir> <stats_output_dir> <bam_files>
            nohup refinem scaffold_stats -c 16 --genome_ext fa ${mypath}/result/spades/spades_${i}/scaffolds_filter_2000.fasta ${mypath}/result/filter_bin/${i} ${mypath}/result/bin_refinem/scaffold_stats/${i} ${mypath}/result/BBmap/${i}.sorted.mapped.bam
        done

#outliers and filter_bins

        mypath=/media/dell/data/huangye/data/meta_all
        mkdir -p ${mypath}/result/bin_refinem/outliers ${mypath}/result/bin_refinem/filter_bins
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do
            #refinem outliers <stats_output_dir>/scaffold_stats.tsv <outlier_output_dir>
            nohup refinem outliers ${mypath}/result/bin_refinem/${i}/scaffold_stats.tsv ${mypath}/result/bin_refinem/outliers/${i}
            #refinem filter_bins <bin_dir> <outlier_output_dir>/outliers.tsv <filtered_output_dir>
            nohup refinem filter_bins --genome_ext fa ${mypath}/result/filter_bin/${i} ${mypath}/result/bin_refinem/outliers/${i}/outliers.tsv ${mypath}/result/bin_refinem/filter_bins/${i}
        done

#4.7 conda activate metawrap1.3.0

        mkdir -p ${mypath}/result/checkm_refinem
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do
            nohup checkm lineage_wf -f ${mypath}/result/checkm_refinem/${i}/bin_checkm.txt -t 16 -x fa ${mypath}/result/bin_refinem/filter_bins/${i} ${mypath}/result/checkm_refinem/${i}
            nohup checkm qa ${mypath}/result/checkm_refinem/${i}/lineage.ms ${mypath}/result/checkm_refinem/${i}
        done

                  
#########################################4.8 bins gene annotations  
                  
#rename

#1. organize genomes and renaming for easier handing

#1.1. copy genomes into a folder named fna

        mkdir -p ${mypath}/bin_results/fna
        for i in `tail -n +2 result/design.txt|cut -f 1`
        do 
            for bin in `ls result/filter_bin/${i}/*.fa | cut -d '/' -f 4`
                do cp ${mypath}/result/bin_refinem/filter_bins/${i}/*fa ${mypath}/bin_results/fna
            done
        done

#change file names

        cd bin_results/fna
        paste <(ls *fa) <(ls *fa | cut -f 1,2 -d '.' | sed 's/\./\_/') | sed 's#^#mv #' > ../script/rename.sh
        sh ../script/rename.sh

#add .fa 

        rename 's/$/.fa/' *

#1.2. rename 

#notice does not need to be done if you have downloaded the genomes, but exemplifies how the genomes were initially named)

        cd fna
        mkdir renamed
        for i in *fa; do awk '/>/{sub(">","&"FILENAME"-");sub(/\.fa/,x)}1' $i > renamed/$i; done
        cd ..

#2. Generate list of files for looping

#Therefore adhere to general rules for renaming:q

#short name, less than 20 characters

#avoid dots

        cd fna/
        ls *fa > Files.txt
        cut -f1 -d "." Files.txt > ../Files.txt
        cd ../..

#protein prediction

#3. call proteins using prokka


#make folder

        mkdir prokka

#run prokka in a loop

        for sample in `cat Files.txt`; do prokka fna/$sample* --outdir prokka/$sample --prefix $sample --locustag $sample --kingdom Archaea --addgenes --force --increment 10 --compliant --centre UU --cpus 20 --norrna --notrna ; done

#4. concatenate files for further analyses

#concatenate

        cat prokka/*/*faa > prokka/All_Bins.faa

#count number of proteins

        grep -c ">" prokka/All_Bins.faa

#clean-up header to remove the annotation

        cut -f1 -d " " prokka/All_Bins.faa > prokka/All_Bins_clean.faa


#4.9 reassembly 

        mkdir -p result/bin_reannotate
        for i in `tail -n +2 result/design.txt|cut -f 1`
        do 
            
            for bin in `ls result/filter_bin/${i}/*.fa | cut -d '/' -f 4`
            do
                    nohup prokka result/bin_reassembly/${i}/${bin} --quiet --rnammer --cpus 48 --outdir result/bin_reannotate/${i}/${bin} --prefix ${bin}
            done
        done

#4.10 eegNOG

###  COG/eggNOG/KEGG

#eggnog4.5

        conda deactivate
        conda activate eggnog

        for i in `tail -n +2 result/design.txt|cut -f 1`
        do 
            
            for bin in `ls result/filter_bin/${i}/*.fa | cut -d '/' -f 4`
            do
                    mkdir -p result/bin_eggnog/${i}/${bin}
                    emapper.py -m diamond --no_annot --no_file_comments --data_dir /share/data1/liushuangjiang/database/eggNOG_4.5_old --cpu 16 -i result/bin_annotate/${i}/${bin}/${bin}.faa -o result/bin_eggnog/${i}/${bin}/${bin}.protein --override

                    emapper.py --annotate_hits_table result/bin_eggnog/${i}/${bin}/${bin}.protein.emapper.seed_orthologs --no_file_comments \
                -o result/bin_eggnog/${i}/${bin}/${bin}.output --cpu 16 --data_dir /share/data1/liushuangjiang/database/eggNOG_4.5_old/ --override
            done
        done


##eggnog5.0

        mkdir -p result/bin_eggnog_5.0
        for i in `tail -n +2 result/design.txt|cut -f 1`
        do 
            
            for bin in `ls result/filter_bin/${i}/*.fa | cut -d '/' -f 4`
            do
                    mkdir -p result/bin_eggnog_5.0/${i}/${bin}
                    emapper.py -m diamond --no_annot --no_file_comments --data_dir /media/dell/data/database/eggNOG5.0 --cpu 16 -i result/bin_annotate/${i}/${bin}/${bin}.faa -o result/bin_eggnog_5.0/${i}/${bin}/${bin}.protein --override

                    emapper.py --annotate_hits_table result/bin_eggnog_5.0/${i}/${bin}/${bin}.protein.emapper.seed_orthologs --no_file_comments -o result/bin_eggnog_5.0/${i}/${bin}/${bin}.output --cpu 16 --data_dir /media/dell/data/database/eggNOG5.0 --override
            done
        done


#eggnog5.0

#'1 i Query_name\teggNOG_ortholog\tEvalue\tScore\tTaxonomic\tProtein\tGO\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tTax_scope\tEggNOG_OGs\tBestOG\tCOG\tDescription'

        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
            for bin in `ls result/filter_bin/${i}/*.fa | cut -d '/' -f 4`
            do
                sed  '1 i Query_name\teggNOG_ortholog\tEvalue\tScore\tTaxonomic\tProtein\tGO\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tTax_scope\tEggNOG_OGs\tBestOG\tCOG\tDescription' result/bin_eggnog_5.0/${i}/${bin}/${bin}.output.emapper.annotations > result/bin_eggnog_pathway_5.0/output/${bin}_output

                cut -f 1,21 result/bin_eggnog_pathway_5.0/output/${bin}_output | cut -f 1 -d ',' | grep -v -P '\t$' > result/bin_eggnog_pathway_5.0/COG/${bin}_1COG.list
                cut -f 1,9  result/bin_eggnog_pathway_5.0/output/${bin}_output|cut -f 1 -d ','|grep -v -P '\t$' > result/bin_eggnog_pathway_5.0/KEGG/${bin}_2ko.list
                cut -f 1,10 result/bin_eggnog_pathway_5.0/output/${bin}_output|cut -f 1 -d ','|grep -v -P '\t$' > result/bin_eggnog_pathway_5.0/KEGG/${bin}_2ko_pathway.list
                cut -f 1,11 result/bin_eggnog_pathway_5.0/output/${bin}_output|cut -f 1 -d ','|grep -v -P '\t$' > result/bin_eggnog_pathway_5.0/KEGG/${bin}_2ko_module.list
                cut -f 1,16 result/bin_eggnog_pathway_5.0/output/${bin}_output|cut -f 1 -d ','|grep -v -P '\t$' > result/bin_eggnog_pathway_5.0/CAZy/${bin}_CAZy.list
                cut -f 1,6,22 result/bin_eggnog_pathway_5.0/output/${bin}_output > result/bin_eggnog_pathway_5.0/${bin}_protein_description.list
            done
        done

### cd-hit create non-rebundant genesets

        mkdir -p ${mypath}/result/NR
        nohup cd-hit-est -i ${mypath}/result/prokka/mg.ffn -o ${mypath}/result/NR/mg.ffn.nr \
          -aS 0.9 -c 0.95 -G 0 -M 0 -T 9 -n 8 -d 0 -g 1
        #ln temp/NR/mg.ffn.nr result/NR/nucleotide.fa
        # convert dna sequences to protein sequences 
        transeq -sequence ${mypath}/result/NR/mg.ffn.nr -outseq ${mypath}/result/NR/protein.fa
        sed -i 's/_1 / /' ${mypath}/result/NR/protein.fa
        ###  salmon
        conda deactivate
        conda activate salmon
        mkdir -p result/salmon
        #nohup salmon index -t result/NR/mg.ffn.nr -p 9 --type quasi -k 31 -i result/salmon/index 
        nohup salmon index -t result/NR/mg.ffn.nr -p 9 -k 31 -i result/salmon/index
        for i in  `tail -n +2 ${mypath}/result/design.txt | cut -f 1`
        do 
                nohup salmon quant -i ${mypath}/result/salmon/index -l A -p 3 --meta --validateMappings -1 ${mypath}/data_cut/${i}_1.fastq -2 ${mypath}/data_cut/${i}_2.fastq -o ${mypath}/result/salmon/${i}.quant
        done
        nohup salmon quantmerge --quants ${mypath}/result/salmon/*.quant -o ${mypath}/result/salmon/gene.TPM
        nohup salmon quantmerge --quants ${mypath}/result/salmon/*.quant --column NumReads -o ${mypath}/result/salmon/gene.count
        sed -i '1 s/.quant//g' ${mypath}/result/salmon/gene.*

##################################################    hmmer   pfam

        conda activate metawrap1.3.0
        nohup hmmscan -o ${mypath}/result/hmmer_gene/hmm.out.txt --tblout ${mypath}/result/hmmer_gene/hmm.out.tbl --noali -E 1e-5 /media/dell/data/database/database_old_machine/pfam/Pfam-A.hmm ${mypath}/result/NR/protein.fa &

##################################################     diamond AMD_nr 

        mypath=/media/dell/data/huangye/data/meta_all
        #
        diamond makedb --in /media/dell/data/huangye/database/AMD_nr_database/nr_AMD/nr.AMD.fa -d /media/dell/data/huangye/database/AMD_nr_database/nr_AMD.diamond
        #
        nohup diamond blastp --db /media/dell/data/huangye/database/AMD_nr_database/nr_AMD.dmnd --query ${mypath}/result/NR/protein.fa --out ${mypath}/result/NR/nr.tab --outfmt 6 --threads 16 --sensitive --max-target-seqs 20 --block-size 16 --evalue 1e-5 --id 30  --tmpdir /dev/shm --index-chunks 1 &

###################################   barrnap  extract 16 RNA #####

        conda activate barrnap
        for i in `tail -n +2 result/design.txt| cut -f 1 `
        do 
            barrnap --kingdom arc ${mypath}/result/spades/spades_${i}/scaffolds.fasta --threads 16 --outseq ${mypath}/result/16S_seq/${i}_arc.fa
            barrnap --kingdom bac ${mypath}/result/spades/spades_${i}/scaffolds.fasta --threads 16 --outseq ${mypath}/result/16S_seq/${i}_bac.fa
        done

###################################  KrakenTools

        for i in  `tail -n +2 result/design.txt | cut -f 1`
        do 
            /media/dell/data/huangye/script/KrakenTools-master/extract_kraken_reads.py -k ${mypath}/result/kraken2/${i}_output -o ${mypath}/result/kraken_exclude/${i}_exclude_bac_1.fastq -o2 ${mypath}/result/kraken_exclude/${i}_exclude_bac_2.fastq -t 179 --exclude --include-children -1 ${mypath}/data_cut/${i}_1.fastq -2 ${mypath}/data_cut/${i}_2.fastq
        done

#################################### ScycDB

        for i in `tail -n +2 design.txt | cut -f 1`;do cp prokka/${i}_prokka/*.ffn prokka/${i}.ffn;done
        mv prokka/*.ffn SCycDB/seq/
        #clean-up header to remove the annotation
        cut -f1 -d " " filter_prokka/All_Bins.faa > filter_prokka/All_Bins_clean.faa
        for i in *fa; do awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' $i > ../for_ycDB/$i; done
        cd /data/huangye/software/SCycDB
        rename 's/.ffn/.fa/' *
        #gene count
         for i in `ls *`;do cat ${i} | grep ">" |wc ;done > ../count.txt

        #sampleinfo.txt 
        1W1  gene count
        #create diamond database
        diamond makedb --in SCycDB_2020Mar.faa -d SCycDB_2020Mar
        #alignment
        perl SCycDB_FunctionProfiler.PL -d for_ycDB -m diamond -f fa -s nucl -si sampleinfo.txt -rs 674734 -o SCycDB
        using SCycDB_TaxonomyProfiler.PL:
        conda activate kraken2
        perl SCycDB_TaxonomyProfiler.PL -d for_ycDB -m diamond -f fa -s nucl -si sampleinfo.txt -rs 674734 

########################################  METABOLIC

        cd /data/huangye/software/METABOLIC
        conda activate METABOLIC
        #METABOLIC-C
        rename 's/.fa/.fasta/' *
        /usr/bin/perl METABOLIC-C.pl -in-gn AMD/result/filter_bin_all -r genome_info.txt -o ./AMD_meta_bin


######################################## Fegenie

        #Fegenie
        conda activate fegenie
        FeGenie.py -h
        for i in `cat bin.list1 | cut -f 1`;do cut -f 1 -d " " faa/${i}.faa > faa/${i}.fa;done
        FeGenie.py -bin_dir /data/huangye/data/MAGs_all/Bins/ffn/ -bin_ext fa -t 16 


######################################## KEGG-decoder

        cd /data/software/kofam_scan-1.3.0
        ./exec_annotation bin_result/all_genome/All_Genomes_clean.faa -f mapper -E 1e-5 --cpu 30 -o hy_MAGs_data/meta_all_kofam_ko.txt
        conda activate KEGGdecoder
        KEGG-decoder --input hy_MAGs_data/meta_all_kofam_ko.txt --output hy_MAGs_data/meta_all_kofam_ko.list --vizoption static


