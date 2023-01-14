
cd /home/yxtu/Data/CRC_Apc/
# 超链接
sed '1d' "/home/yxtu/Data/CRC_Apc/sampleInfo.csv" |awk -F "," '{print $1"\t"$3}'|xargs -n2 bash -c 'ln -s /home/kcyu/project/Apc_metagenomic/$0/2*/*_1.fq.gz $1_1.fq.gz'
sed '1d' "/home/yxtu/Data/CRC_Apc/sampleInfo.csv" |awk -F "," '{print $1"\t"$3}'|xargs -n2 bash -c 'ln -s /home/kcyu/project/Apc_metagenomic/$0/2*/*_2.fq.gz $1_2.fq.gz'

mkdir -p /home/yxtu/Data/CRC_Apc/kneaddata/output

for i in $(ls /home/yxtu/Data/CRC_Apc/*_1.fq.gz); do j=`basename $i`;
echo "kneaddata -i ${i} -i ${i%%_1.fq.gz}_2.fq.gz -o /home/yxtu/Data/CRC_Apc/kneaddata/output/ -v \
--reference-db /home/yxtu/Data/CRC_Apc/kneaddata/mouse_genome/  \
--trimmomatic /home/yxtu/software/Trimmomatic-0.38/  \
--trimmomatic-options 'SLIDINGWINDOW:4:20 MINLEN:50' \
-t 9 --bowtie2-options '--very-sensitive --dovetail' --remove-intermediate-output" > /home/yxtu/Data/CRC_Apc/kneaddata/output/${j%%_1.fq.gz}.pbs;done
for i in $(ls /home/yxtu/Data/CRC_Apc/kneaddata/output/*.pbs); do qsub -q batch -V -l nodes=2:ppn=4 $i;done


for i in $(ls /home/yxtu/Data/CRC_Apc/kneaddata/output/*_paired_1.fastq); do j=`basename $i`;
echo "cat  ${i} ${i%%_paired_1.fastq}_paired_2.fastq ${i%%_paired_1.fastq}_unmatched_1.fastq ${i%%_paired_1.fastq}_unmatched_2.fastq  | awk '{if(NR%4==1) print \"@\"NR; else print \$0}'> /home/yxtu/Data/CRC_Apc/kneaddata/cat_reads/${j%%_paired_1.fastq}.fastq" >/home/yxtu/Data/CRC_Apc/kneaddata/${j%%_paired_1.fastq}.pbs;done
for i in $(ls  /home/yxtu/Data/CRC_Apc/kneaddata/cat_reads/*.pbs); do qsub -q batch -V -l nodes=2:ppn=4 $i;done

for i in $(ls   /home/yxtu/Data/CRC_Apc/kneaddata/cat_reads/*.fastq); do j=`basename $i`;
echo "metaphlan ${i} --input_type fastq -o /home/yxtu/Data/CRC_Apc/metaphlan/${j%%.fastq}_metaphlan.txt" > /home/yxtu/Data/CRC_Apc/metaphlan/${j%%.fastq}_metaphlan.pbs;done
for i in $(ls /home/yxtu/Data/CRC_Apc/metaphlan/*_metaphlan.pbs); do qsub -q batch -V -l nodes=2:ppn=4 $i;done

merge_metaphlan_tables.py /home/yxtu/Data/CRC_Apc/metaphlan/*.txt  >/home/yxtu/Data/CRC_Apc/metaphlan/merged_abundance_table.txt
cd /home/yxtu/Data/CRC_Apc/metaphlan/
grep -E '(s__)|(clade_name)' merged_abundance_table.txt |grep -v 't__'|sed 's/^.*s__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' > merged_abundance_table_species.txt
grep -E '(g__)|(clade_name)' merged_abundance_table.txt |grep -v 's__'|sed 's/^.*g__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' > merged_abundance_table_genus.txt
grep -E '(f__)|(clade_name)' merged_abundance_table.txt |grep -v 'g__'|sed 's/^.*f__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' > merged_abundance_table_family.txt
grep -E '(o__)|(clade_name)' merged_abundance_table.txt |grep -v 'f__'|sed 's/^.*o__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' > merged_abundance_table_order.txt
grep -E '(c__)|(clade_name)' merged_abundance_table.txt |grep -v 'o__'|sed 's/^.*c__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' > merged_abundance_table_class.txt
grep -E '(p__)|(clade_name)' merged_abundance_table.txt |grep -v 'c__'|sed 's/^.*p__//g'|awk '{$2=null;print}'|sed 's/\ \ /\ /g'|sed 's/\ /\t/g' > merged_abundance_table_phylum.txt

# stain-level
mkdir /home/yxtu/Data/CRC_Apc/strain/
cd /home/yxtu/Data/CRC_Apc/strain/
mkdir -p sams
mkdir -p bowtie2
mkdir -p profiles

for i in $(ls /home/yxtu/Data/CRC_Apc/kneaddata/cat_reads/*.fastq); do j=`basename $i`;
echo "metaphlan ${i} --input_type fastq -s /home/yxtu/Data/CRC_Apc/strain/sams/${j%%.fastq}.sam.bz2 --bowtie2out /home/yxtu/Data/CRC_Apc/strain/bowtie2/${j%%.fastq}.bowtie2.bz2 -o /home/yxtu/Data/CRC_Apc/strain/profiles/${j%%.fastq}_profiled.tsv" > /home/yxtu/Data/CRC_Apc/strain/step_001_${j%%.fastq}.pbs
done
for i in $(ls /home/yxtu/Data/CRC_Apc/strain/step_001*.pbs); do qsub -q batch -V -l nodes=2:ppn=4 $i;done

mkdir -p /home/yxtu/Data/CRC_Apc/strain/consensus_markers
sample2markers.py -i /home/yxtu/Data/CRC_Apc/strain/sams/*.sam.bz2 -o /home/yxtu/Data/CRC_Apc/strain/consensus_markers -n 8


mkdir -p /home/yxtu/Data/CRC_Apc/strain/db_markers
extract_markers.py -d /home/yxtu/miniconda3/envs/metaphlan/lib/python3.7/site-packages/metaphlan/utils/../metaphlan_databases/mpa_v30_CHOCOPhlAn_201901.pkl -c s__Akkermansia_muciniphila -o /home/yxtu/Data/CRC_Apc/strain/db_markers/

mkdir -p /home/yxtu/Data/CRC_Apc/strain/output
strainphlan -s /home/yxtu/Data/CRC_Apc/strain/consensus_markers/*.pkl -m /home/yxtu/Data/CRC_Apc/strain/db_markers/s__Akkermansia_muciniphila.fna -o /home/yxtu/Data/CRC_Apc/strain/output -n 8 -c s__Akkermansia_muciniphila --mutation_rates


# add_metadata_tree.py -t "/home/yxtu/Data/CRC_Apc/strain/output/RAxML_bestTree.s__Akkermansia_muciniphila.StrainPhlAn3.tre" -f /home/yxtu/Data/CRC_Apc/strain/output/metadata.txt 
conda activate py27
miniconda3/envs/biobakry/bin/plot_tree_graphlan.py -t "/home/yxtu/Data/CRC_Apc/strain/output/RAxML_bestTree.s__Akkermansia_muciniphila.StrainPhlAn3.tre.metadata" 




# humann3
for i in $(ls /home/yxtu/Data/CRC_Apc/kneaddata/cat_reads/*.fastq);do j=`basename $i`;
echo "humann --threads 9 --input ${i}  --output /home/yxtu/Data/CRC_Apc/humann2_out/${j%%.fastq}/" > /home/yxtu/Data/CRC_Apc/humann2_out/${j%%.fastq}.pbs;done
for i in $(ls /home/yxtu/Data/CRC_Apc/humann2_out/*.pbs); do qsub -q batch -V -l nodes=2:ppn=6 $i;done

mkdir "/home/yxtu/Data/CRC_Apc/humann2_out/genefamily/"
mkdir "/home/yxtu/Data/CRC_Apc/humann2_out/pathway/"
mv /home/yxtu/Data/CRC_Apc/humann2_out/*_kneaddata/*_kneaddata_pathabundance.tsv "/home/yxtu/Data/CRC_Apc/humann2_out/pathway/"
mv /home/yxtu/Data/CRC_Apc/humann2_out/*_kneaddata/*_kneaddata_genefamilies.tsv "/home/yxtu/Data/CRC_Apc/humann2_out/genefamily/"

# 合并文件
mkdir /home/yxtu/Data/CRC_Apc/humann2_out/result/
humann_join_tables --input /home/yxtu/Data/CRC_Apc/humann2_out/genefamily/ --output /home/yxtu/Data/CRC_Apc/humann2_out/result/humann2_genefamilies.tsv 
humann_join_tables --input /home/yxtu/Data/CRC_Apc/humann2_out/pathway/ --output /home/yxtu/Data/CRC_Apc/humann2_out/result/humann2_pathways.tsv 
# 基因注释
for i in {rxn,go,ko,level4ec,pfam,eggnog};
do humann_regroup_table --input /home/yxtu/Data/CRC_Apc/humann2_out/result/humann2_genefamilies.tsv  --groups uniref90_${i} --output /home/yxtu/Data/CRC_Apc/humann2_out/result/humann2_anno_${i}.tsv;
done

#标准化归一化
for i in $(ls /home/yxtu/Data/CRC_Apc/humann2_out/result/*.tsv);
do humann_renorm_table --input ${i} --output ${i%%.tsv}_relab.tsv --units relab --special n;done

for i in $(ls /home/yxtu/Data/CRC_Apc/humann2_out/result/*.tsv);
do humann_renorm_table --input ${i} --output ${i%%.tsv}_cpm.tsv --units cpm --special n;done