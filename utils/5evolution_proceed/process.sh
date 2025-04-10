基因课全部课程
通过网盘分享的文件：【18天大课】基因组与比较基因组
链接: https://pan.baidu.com/s/1tZcBGRyOoCrwLlD217Vq_g?pwd=39sz 提取码: 39sz 
--来自百度网盘超级会员v6的分享
#复制备份gff
cp Oryza_indica.ASM465v1.60.gff3 Oryza_indica.gff3.tmp
cp Oryza_sativa.IRGSP-1.0.60.gff3 Oryza_sativa.gff3.tmp
cp Sorghum_bicolor.Sorghum_bicolor_NCBIv3.60.gff3 Sorghum_bicolor.gff3.tmp
cp Triticum_aestivum.IWGSC.60.gff3 Triticum_aestivum.gff3.tmp
cp Zea_mays.Zm-B73-REFERENCE-NAM-5.0.60.gff3 Zea_mays.gff3.tmp

#替换gene:
sudo sed -i 's/=gene:/=/g' Oryza_indica.gff3.tmp
sudo sed -i 's/=transcript:/=/g' Oryza_indica.gff3.tmp
sudo sed -i 's/=CDS:/=/g' Oryza_indica.gff3.tmp

sudo sed -i 's/=gene:/=/g' Oryza_sativa.gff3.tmp
sudo sed -i 's/=transcript:/=/g' Oryza_sativa.gff3.tmp
sudo sed -i 's/=CDS:/=/g' Oryza_sativa.gff3.tmp

sudo sed -i 's/=gene:/=/g' Sorghum_bicolor.gff3.tmp
sudo sed -i 's/=transcript:/=/g' Sorghum_bicolor.gff3.tmp
sudo sed -i 's/=CDS:/=/g' Sorghum_bicolor.gff3.tmp

sudo sed -i 's/=gene:/=/g' Triticum_aestivum.gff3.tmp
sudo sed -i 's/=transcript:/=/g' Triticum_aestivum.gff3.tmp
sudo sed -i 's/=CDS:/=/g' Triticum_aestivum.gff3.tmp

sudo sed -i 's/=gene:/=/g' Zea_mays.gff3.tmp
sudo sed -i 's/=transcript:/=/g' Zea_mays.gff3.tmp
sudo sed -i 's/=CDS:/=/g' Zea_mays.gff3.tmp

#提取最长转录本
perl /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/gff_ensembl_longest.pl Oryza_indica.gff3.tmp Oryza_indica_longest.gene2mrna_id Oryza_indica.gff3
perl /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/gff_ensembl_longest.pl Oryza_sativa.gff3.tmp Oryza_sativa_longest.gene2mrna_id Oryza_sativa.gff3
perl /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/gff_ensembl_longest.pl Sorghum_bicolor.gff3.tmp Sorghum_bicolor_longest.gene2mrna_id Sorghum_bicolor.gff3
perl /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/gff_ensembl_longest.pl Triticum_aestivum.gff3.tmp Triticum_aestivum_longest.gene2mrna_id Triticum_aestivum.gff3
perl /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/gff_ensembl_longest.pl Zea_mays.gff3.tmp Zea_mays_longest.gene2mrna_id Zea_mays.gff3

#提取mRNAid
awk '{print $2}' Oryza_indica_longest.gene2mrna_id > Oryza_indica.mrna_id
awk '{print $2}' Oryza_sativa_longest.gene2mrna_id > Oryza_sativa.mrna_id
awk '{print $2}' Sorghum_bicolor_longest.gene2mrna_id > Sorghum_bicolor.mrna_id
awk '{print $2}' Triticum_aestivum_longest.gene2mrna_id > Triticum_aestivum.mrna_id
awk '{print $2}' Zea_mays_longest.gene2mrna_id > Zea_mays.mrna_id

#提取cds和pep
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.ASM465v1.cds.all.fa Oryza_indica.mrna_id  >  Oryza_indica.cds.fa
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.ASM465v1.pep.all.fa Oryza_indica_pep.mrna_id  >  Oryza_indica.pep.fa

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.IRGSP-1.0.cds.all.fa Oryza_sativa.mrna_id  >  Oryza_sativa.cds.fa
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.IRGSP-1.0.pep.all.fa Oryza_sativa.mrna_id  >  Oryza_sativa.pep.fa

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.cds.all.fa Sorghum_bicolor.mrna_id  >  Sorghum_bicolor.cds.fa
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa Sorghum_bicolor.mrna_id  >  Sorghum_bicolor.pep.fa

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.IWGSC.cds.all.fa Triticum_aestivum.mrna_id  >  Triticum_aestivum.cds.fa
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.IWGSC.pep.all.fa Triticum_aestivum.mrna_id  >  Triticum_aestivum.pep.fa

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cds.all.fa Zea_mays.mrna_id  >  Zea_mays.cds.fa
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif seqtk subseq /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.pep.all.fa Zea_mays_pep.mrna_id  >  Zea_mays.pep.fa

#基因家族聚类
#https://orthovenn3.bioinfotoolkits.net/start/db 在线处理
#https://orthovenn3.bioinfotoolkits.net/result/5b369d6b01c54be0a2f55b79c26a81e7/orthologous  在线处理

#数据链接
cd /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/2GeneFamilyCluster/data
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.pep.fa indica.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.pep.fa japonica.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.pep.fa sorghum.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.pep.fa wheat.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.pep.fa maize.fasta


#基因家族聚类
cd /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/2GeneFamilyCluster
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif orthofinder  -f data \
	-S diamond \
	-M msa \
	-T fasttree \
	-t 20

#韦恩分析
cp /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/2GeneFamilyCluster/data/OrthoFinder/Results_Nov01/Orthogroups/Orthogroups.GeneCount.tsv ./
dos2unix Orthogroups.GeneCount.tsv
awk '{ if(NR==1){ for(i=2;i<NF;i++ ){printf $i"\t"} }else{  for(i=2;i<NF;i++){ if($i>0){printf $1 } ; printf "\t" }}; printf "\n"}' Orthogroups.GeneCount.tsv|sed 's/\t$//' > Orthogroups.GeneCount.venn


#系统发育树构建

cp /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/2GeneFamilyCluster/data/OrthoFinder/Results_Nov01/MultipleSequenceAlignments/SpeciesTreeAlignment.fa /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/3PhylogeneticTree
cp /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/02.PhylogeneticTree/test2_Orthofinder/run_tree.sh /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/3PhylogeneticTree

#修剪比对结果
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  trimal \
    -in  SpeciesTreeAlignment.fa  \
    -out SpeciesTreeAlignment_trim.fa \
    -fasta \
    -gt 0.6  \
    -cons 60

#构建物种树
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif raxmlHPC-PTHREADS  \
    -T 20 \
    -m PROTGAMMAJTT \
    -f a \
    -p 123 -x 123 -# 100 \
    -n out \
    -s  SpeciesTreeAlignment_trim.fa   1>tree.log 2>tree.err

#整理Astral输入
orthDir=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/2GeneFamilyCluster/data/OrthoFinder/Results_Nov01/
cat  $orthDir/Orthogroups/Orthogroups_SingleCopyOrthologues.txt |while read aa; do cat  $orthDir/Gene_Trees/$aa\_tree.txt |awk '{print $0 }'   ;done > SingleCopy.trees

sed -r  's/([(,]{1}[A-Za-z]+)_[^:]+/\1/g' SingleCopy.trees > Astral_input.trees

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif java -jar   /opt/ASTRAL/astral.jar  -i  Astral_input.trees  -o Astral_output.tree 2>out.log

#在线展示进化树的方式
#MEGA Figtree  在线网站 Evolview iTOL

#分歧时间估计  分化时间估计
#mega timetree   http://www.timetree.org/home


##=================基于mcmc氨基酸序列====================
orthoDir=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/2GeneFamilyCluster/data/OrthoFinder/Results_Nov01/
## 拷贝进化树
cp $orthoDir/Species_Tree/SpeciesTree_rooted.txt ./
## 拷贝多序列比对结果
cp $orthoDir/MultipleSequenceAlignments/SpeciesTreeAlignment.fa ./

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  trimal  -in SpeciesTreeAlignment.fa -out supergene.phy -phylip_paml

sed 's/:[^,)(]\+//g' SpeciesTree_rooted.txt|sed 's/)1/)/g' > input.tree
vi input.tree
 5  1


singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  mcmctree mcmctree.ctl
cp tmp0001.ctl codeml.ctl
echo "clock = 1" >> codeml.ctl
sed -i 's/\(getSE.*=\s*\d/getSE = 0/' codeml.ctl
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  codeml  codeml.ctl

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif mcmctree  mcmctree.ctl

mv out.BV in.BV

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif mcmctree  mcmctree.ctl 

#=============================基于4D位点=================================

orthoDir=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/2GeneFamilyCluster/data/OrthoFinder/Results_Nov01/

cp   $orthoDir/Orthogroups/Orthogroups.tsv ./
cp   $orthoDir/Orthogroups/Orthogroups_SingleCopyOrthologues.txt ./

dos2unix  Orthogroups.tsv

# 生成单拷贝基因家族列表文件

## 提取单拷贝基因家族成员方法1  
awk  '{if( NR==FNR ){ A[$1]=1 } else{ if($1 in A || FNR==1){print $0}} }'  Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > single_copy.txt

## 提取单拷贝基因家族成员方法2
awk  '{if(NR==1){f=NF}; if( $0 !~ /,/ && NF==f ){print $0}}' Orthogroups.tsv > single_copy.txt


## 对基因ID添加物种前缀
awk  '{if(NR==1){ n=split($0,A, "\t")} else { for(i=2; i<= n; i++){  printf A[i]"_"$i"\t" }; printf "\n" } }' single_copy.txt |sed 's/\s\+$//' > single_copy.txt.change


# 生成cds水平多序列比对结果

## 生成ParaAT线程数文件
echo "6" > proc.txt

ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.cds.fa ./Oryza_indica.cds.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.pep.fa ./Oryza_indica.pep.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.cds.fa ./Oryza_sativa.cds.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.pep.fa ./Oryza_sativa.pep.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.cds.fa ./Sorghum_bicolor.cds.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.pep.fa ./Sorghum_bicolor.pep.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.cds.fa ./Triticum_aestivum.cds.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.pep.fa ./Triticum_aestivum.pep.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.cds.fa ./Zea_mays.cds.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.pep.fa ./Zea_mays.pep.fa

sed 's/-TA/-PA/g' Oryza_indica.cds.fa > Oryza_indica_modified.cds.fa
sed 's/_T00/_P00/g' Zea_mays.cds.fa > Zea_mays_modified.cds.fa

awk '/^>/ {print $1; next} {print}' Oryza_indica.pep.fa > formatted_indica.pep.fa
awk '/^>/ {print $1; next} {print}' Oryza_indica_modified.cds.fa > formatted_indica.cds.fa
awk '/^>/ {print $1; next} {print}' Oryza_sativa.pep.fa > formatted_japonica.pep.fa
awk '/^>/ {print $1; next} {print}' Oryza_sativa.cds.fa > formatted_japonica.cds.fa
awk '/^>/ {print $1; next} {print}' Sorghum_bicolor.pep.fa > formatted_sorghum.pep.fa
awk '/^>/ {print $1; next} {print}' Sorghum_bicolor.cds.fa > formatted_sorghum.cds.fa
awk '/^>/ {print $1; next} {print}' Triticum_aestivum.pep.fa > formatted_wheat.pep.fa
awk '/^>/ {print $1; next} {print}' Triticum_aestivum.cds.fa > formatted_wheat.cds.fa
awk '/^>/ {print $1; next} {print}' Zea_mays.pep.fa > formatted_maize.pep.fa
awk '/^>/ {print $1; next} {print}' Zea_mays_modified.cds.fa > formatted_maize.cds.fa


sed 's/^>/>indica_/' formatted_indica.cds.fa  > Oryza_indica.cds.fa.change
sed 's/^>/>indica_/' formatted_indica.pep.fa  > Oryza_indica.pep.fa.change
sed 's/^>/>japonica_/' formatted_japonica.cds.fa  > Oryza_sativa.cds.fa.change
sed 's/^>/>japonica_/' formatted_japonica.pep.fa  > Oryza_sativa.pep.fa.change
sed 's/^>/>sorghum_/' formatted_sorghum.cds.fa  > Sorghum_bicolor.cds.fa.change
sed 's/^>/>sorghum_/' formatted_sorghum.pep.fa  > Sorghum_bicolor.pep.fa.change
sed 's/^>/>wheat_/' formatted_wheat.cds.fa  > Triticum_aestivum.cds.fa.change
sed 's/^>/>wheat_/' formatted_wheat.pep.fa  > Triticum_aestivum.pep.fa.change
sed 's/^>/>maize_/' formatted_maize.cds.fa  > Zea_mays.cds.fa.change
sed 's/^>/>maize_/' formatted_maize.pep.fa  > Zea_mays.pep.fa.change

cat *.cds.fa.change > all.cds.fa
cat *.pep.fa.change > all.pep.fa


## 生成cds水平比对结果
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif ParaAT_mdf.pl -h single_copy.txt.change -a all.pep.fa -n all.cds.fa   -o Para_out -p  proc.txt  

# 合并单个家族成supergene phylip格式文件

## 去掉基因ID，仅保留物种简写
ls Para_out/*.cds_aln.fasta |while read aa;do awk -F "_" '{print $1}' $aa > $aa.change;done

## 对单个家族进行连接
seqkit concat Para_out/*.cds_aln.fasta.change  >  single_copy.cds_msa.fasta

## fasta格式转phylip格式
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  trimal  -in single_copy.cds_msa.fasta -out single_copy.cds_msa.phy -phylip_paml 


# 提取4d位点
perl  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/extract_4d_phy.pl  single_copy.cds_msa.phy  single_copy.cds_msa.4d.phy

# 进行分化时间计算
singularity exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif   mcmctree mcmctree.ctl

# 以下为手动生成单拷贝基因家族supergene蛋白序列方法 
ls Para_out/*.pep_aln  |while read aa;do awk -F "_" '{print $1}' $aa > $aa.change;done

seqkit concat Para_out/*.pep_aln.change  >  single_copy.pep_msa.fasta

singularity  exec  ../../software/PhyloTools.sif  trimal  -in single_copy.pep_msa.fasta -out single_copy.pep_msa.phy -phylip_paml 

#==============================基因家族收缩与扩张=================================
dos2unix Orthogroups.GeneCount.tsv
sed 's/[a-zA-Z0-9]\+$//' Orthogroups.GeneCount.tsv | awk '{print $1"\t"$0}'   > input.tab

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  cafe5  --infile input.tab --tree input.tree --output_prefix  cafe_ortho  --cores 5


#整理绘图饼图
sed -r 's/^([^>]+)<[0-9]+>/\1/' Base_clade_results.txt

grep "TREE " cafe_ortho/Base_asr.tre |wc -l
#19511

sed -r 's/^([^>]+)<[0-9]+>/\1/' Base_clade_results.txt | awk '!/^#/{ p=0.5; s=10; C=19511-$2-$3; if(!/^</){p=-1}; print $1 "\t" p "\t" s "\t" $2 "\t" $3 "\t" C }'

cat cafe_ortho/Base_clade_results.txt|awk '{print $1"\t+"$2"/-"$3"\t0.5\t#ff0000\tnormal\t1\t0"}'

#============================共线性分析===========================================

#准备数据

ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.pep.fa ./Oryza_indica.pep.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.pep.fa ./Oryza_sativa.pep.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.pep.fa ./Sorghum_bicolor.pep.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.pep.fa ./Triticum_aestivum.pep.fa
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.pep.fa ./Zea_mays.pep.fa
sed 's/^>/>indica|/' Oryza_indica.pep.fa  > indica.fasta
sed 's/^>/>japonica|/'  Oryza_sativa.pep.fa > japonica.fasta
sed 's/^>/>sorghum|/' Sorghum_bicolor.pep.fa  > sorghum.fasta
sed 's/^>/>wheat|/'  Triticum_aestivum.pep.fa > wheat.fasta
sed 's/^>/>maize|/' Zea_mays.pep.fa  > maize.fasta

cat  indica.fasta japonica.fasta sorghum.fasta wheat.fasta maize.fasta > all.fasta

makeblastdb -in  all.fasta -dbtype prot

blastp -query  all.fasta -db all.fasta -out blast_out.tab  -outfmt 6 -evalue 1e-10  -num_threads 30

#物种自身与自身
awk '$1~"indica" && $2~"indica" && $1 != $2'  blast_out.tab | sed 's/indica|//g' >indica.blast
awk '$1~"japonica" && $2~"japonica" && $1 != $2'  blast_out.tab | sed 's/japonica|//g' >japonica.blast
awk '$1~"wheat" && $2~"wheat" && $1 != $2'  blast_out.tab | sed 's/wheat|//g' >wheat.blast
awk '$1~"maize" && $2~"maize" && $1 != $2'  blast_out.tab | sed 's/maize|//g' >maize.blast
awk '$1~"sorghum" && $2~"sorghum" && $1 != $2'  blast_out.tab | sed 's/sorghum|//g' >sorghum.blast

#物种与物种之间
awk '$1~"indica" && $2~"japonica"' blast_out.tab | sed 's/indica|//g' | sed 's/japonica|//g' >japonica_indica.blast
awk '$1~"japonica" && $2~"indica"' blast_out.tab | sed 's/indica|//g' | sed 's/japonica|//g' >>japonica_indica.blast 

awk '$1~"indica" && $2~"wheat"' blast_out.tab | sed 's/indica|//g' | sed 's/wheat|//g' >indica_wheat.blast
awk '$1~"wheat" && $2~"indica"' blast_out.tab | sed 's/wheat|//g' | sed 's/indica|//g' >>indica_wheat.blast 

awk '$1~"wheat" && $2~"maize"' blast_out.tab | sed 's/wheat|//g' | sed 's/maize|//g' >wheat_maize.blast
awk '$1~"maize" && $2~"wheat"' blast_out.tab | sed 's/maize|//g' | sed 's/wheat|//g' >>wheat_maize.blast 

awk '$1~"maize" && $2~"sorghum"' blast_out.tab | sed 's/maize|//g' | sed 's/sorghum|//g' >maize_sorghum.blast
awk '$1~"sorghum" && $2~"maize"' blast_out.tab | sed 's/sorghum|//g' | sed 's/maize|//g' >>maize_sorghum.blast 


ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.gff3 ./Oryza_indica.gff3
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.gff3 ./Oryza_sativa.gff3
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.gff3 ./Sorghum_bicolor.gff3
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.gff3 ./Triticum_aestivum.gff3
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.gff3 ./Zea_mays.gff3

less Oryza_indica.gff3|cut -f 1|sort -u|less
grep "^[0-9]" Oryza_indica.gff3 |awk -F '\t|;' '$3=="mRNA"{print $1"\t"$9"\t"$4"\t"$5}' | sed 's/ID=//' >indica.gff
less Oryza_sativa.gff3|cut -f 1|sort -u|less
grep "^[0-9]" Oryza_sativa.gff3 |awk -F '\t|;' '$3=="mRNA"{print $1"\t"$9"\t"$4"\t"$5}' | sed 's/ID=//' >japonica.gff
less Sorghum_bicolor.gff3|cut -f 1|sort -u|less
grep "^[0-9]" Sorghum_bicolor.gff3 |awk -F '\t|;' '$3=="mRNA"{print $1"\t"$9"\t"$4"\t"$5}' | sed 's/ID=//' >sorghum.gff
less Triticum_aestivum.gff3|cut -f 1|sort -u|less
grep "^[0-9]" Triticum_aestivum.gff3 |awk -F '\t|;' '$3=="mRNA"{print $1"\t"$9"\t"$4"\t"$5}' | sed 's/ID=//' >wheat.gff
less Zea_mays.gff3|cut -f 1|sort -u|less
grep "^[0-9]" Zea_mays.gff3 |awk -F '\t|;' '$3=="CDS"{print $1"\t"$9"\t"$4"\t"$5}' | sed 's/ID=//' >maize.gff

cat indica.gff  japonica.gff >indica_japonica.gff

cat japonica.gff indica.gff >japonica_indica.gff 
cat indica.gff  wheat.gff >indica_wheat.gff 
cat wheat.gff  maize.gff >wheat_maize.gff 
cat maize.gff  sorghum.gff >maize_sorghum.gff  

# 共线性分析

## 种内
singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif  MCScanX   indica
singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif  MCScanX   japonica
singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif  MCScanX   wheat
singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif  MCScanX   maize
singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif  MCScanX   sorghum



## 种间
prefix=maize_sorghum #japonica_indica indica_wheat wheat_maize 
singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif  MCScanX   $prefix
# 共线性画图

singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif  java dot_plotter -g $prefix.gff -s $prefix.collinearity -c $prefix.dot.ctl -o $prefix.dot.png

singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif java dual_synteny_plotter -g $prefix.gff -s $prefix.collinearity -c $prefix.dot.ctl -o $prefix.dual_synteny.png

singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif java circle_plotter -g $prefix.gff -s $prefix.collinearity -c $prefix.dot.ctl -o $prefix.circle.png

singularity  exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif java bar_plotter -g $prefix.gff -s $prefix.collinearity -c $prefix.dot.ctl -o $prefix.bar.png


prefix_list=("maize_sorghum" "japonica_indica" "indica_wheat" "wheat_maize")

# 运行循环
for prefix in "${prefix_list[@]}"
do
    echo "Processing $prefix..."

    # MCScanX
    singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif MCScanX $prefix

    # 共线性画图
    singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif java dot_plotter -g $prefix.gff -s $prefix.collinearity -c $prefix.dot.ctl -o $prefix.dot.png

    singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif java dual_synteny_plotter -g $prefix.gff -s $prefix.collinearity -c $prefix.dot.ctl -o $prefix.dual_synteny.png

    singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif java circle_plotter -g $prefix.gff -s $prefix.collinearity -c $prefix.dot.ctl -o $prefix.circle.png

    singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif java bar_plotter -g $prefix.gff -s $prefix.collinearity -c $prefix.dot.ctl -o $prefix.bar.png

    echo "$prefix processing complete."
done

#######=========================用jcvi==============================
## 过滤和格式调整
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.cds.fa  indica.cds
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.cds.fa  japonica.cds
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.cds.fa  sorghum.cds
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.cds.fa  wheat.cds
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.cds.fa  maize.cds

ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.gff3 ./Oryza_indica.gff3
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.gff3 ./Oryza_sativa.gff3
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.gff3 ./Sorghum_bicolor.gff3
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.gff3 ./Triticum_aestivum.gff3
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.gff3 ./Zea_mays.gff3

grep '^[0-9]' Oryza_indica.gff3 > indica.filter.gff3
grep '^[0-9]' Oryza_sativa.gff3 > japonica.filter.gff3
grep '^[0-9]' Sorghum_bicolor.gff3 > sorghum.filter.gff3
grep '^[0-9]' Triticum_aestivum.gff3 > wheat.filter.gff3
grep '^[0-9]' Zea_mays.gff3 > maize.filter.gff3

## 生成 bed 
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.formats.gff bed --type=mRNA indica.filter.gff3 -o indica.bed
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.formats.gff bed --type=mRNA japonica.filter.gff3 -o japonica.bed
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.formats.gff bed --type=mRNA sorghum.filter.gff3 -o sorghum.bed
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.formats.gff bed --type=mRNA wheat.filter.gff3 -o wheat.bed
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.formats.gff bed --type=mRNA maize.filter.gff3 -o maize.bed


# 共线性模块鉴定
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.catalog ortholog japonica indica --no_strip_names
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.catalog ortholog indica wheat --no_strip_names
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.catalog ortholog wheat maize --no_strip_names
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.catalog ortholog maize sorghum --no_strip_names


# 可视化
## 点图
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.graphics.dotplot -o japonica.indica.anchors.pdf  japonica.indica.anchors
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.graphics.dotplot -o indica.wheat.anchors.pdf  indica.wheat.anchors
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.graphics.dotplot -o wheat.maize.anchors.pdf  wheat.maize.anchors
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.graphics.dotplot -o maize.sorghum.anchors.pdf  maize.sorghum.anchors

## 共线性图
### 两个物种
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.synteny screen --minspan=30 --simple japonica.indica.anchors japonica.indica.anchors.new
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.synteny screen --minspan=30 --simple indica.wheat.anchors indica.wheat.anchors.new
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.synteny screen --minspan=30 --simple wheat.maize.anchors wheat.maize.anchors.new
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.compara.synteny screen --minspan=30 --simple maize.sorghum.anchors maize.sorghum.anchors.new


singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.graphics.karyotype seqids2 layout2

### 三个物种
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.graphics.karyotype seqids3 layout3

##5个物种
singularity exec /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/Collinearity.sif python3 -m jcvi.graphics.karyotype seqids5 layout5


#==========================kaks计算=======================
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.cds.fa  indica.cds.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.cds.fa  japonica.cds.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.cds.fa  sorghum.cds.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.cds.fa  wheat.cds.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.cds.fa  maize.cds.fasta

ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_indica.pep.fa  indica.pep.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Oryza_sativa.pep.fa  japonica.pep.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Sorghum_bicolor.pep.fa  sorghum.pep.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Triticum_aestivum.pep.fa  wheat.pep.fasta
ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/1Genome/Zea_mays.pep.fa  maize.pep.fasta

cds1=indica.cds.fasta
cds2=japonica.cds.fasta
cds3=sorghum.cds.fasta
cds4=wheat.cds.fasta
cds5=maize.cds.fasta


## diamond 比对 并筛选直系同源基因对
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif wgd   dmd --eval 1e-10  -o japonica_indica_dmd_out --nostrictcds japonica.cds.fasta indica.cds.fasta
cp japonica_indica_dmd_out/*.rbh ./japonica_indica_homo_pairs.txt

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif wgd   dmd --eval 1e-10  -o indica_wheat_dmd_out --nostrictcds indica.cds.fasta wheat.cds.fasta
cp indica_wheat_dmd_out/*.rbh ./indica_wheat_homo_pairs.txt

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif wgd   dmd --eval 1e-10  -o wheat_maize_dmd_out --nostrictcds wheat.cds.fasta maize.cds.fasta
cp wheat_maize_dmd_out/*.rbh ./wheat_maize_homo_pairs.txt

singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif wgd   dmd --eval 1e-10  -o maize_sorghum_dmd_out --nostrictcds maize.cds.fasta sorghum.cds.fasta
cp maize_sorghum_dmd_out/*.rbh ./maize_sorghum_homo_pairs.txt


#========================================计算============================================
sed -i 's/-PA/-TA/g' indica.pep.fasta
# grep "BGIOSGA012944-TA" indica.pep.fasta

# grep "Os03t0436300-00" japonica.cds.fasta

#---------------------------------japoncia-indica示例----------------------------------
awk '/^>/ {print $1; next} {print}' indica.pep.fasta > formatted_indica.pep.fasta
awk '/^>/ {print $1; next} {print}' indica.cds.fasta > formatted_indica.cds.fasta
awk '/^>/ {print $1; next} {print}' japonica.pep.fasta > formatted_japonica.pep.fasta
awk '/^>/ {print $1; next} {print}' japonica.cds.fasta > formatted_japonica.cds.fasta

cat formatted_japonica.cds.fasta formatted_indica.cds.fasta > input_japonica_indica_cds.fasta
cat formatted_japonica.pep.fasta formatted_indica.pep.fasta > input_japonica_indica_pep.fasta

## 数据量比较大，提取前100对ortholog进行分分析
head -n 10 japonica_indica_homo_pairs.txt > japonica_indica_homo_pairs.txt1

## 蛋白比对转cds
echo "6" > proc.txt
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  ParaAT_mdf.pl -h japonica_indica_homo_pairs.txt  -a input_japonica_indica_pep.fasta -n input_japonica_indica_cds.fasta -p proc.txt  -o  japonica_indica_align_out  -m muscle -f axt

cat japonica_indica_align_out/*.axt > japonica_indica_merge_align.axt

## calculate kaks 
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif KaKs_Calculator -m YN -i japonica_indica_merge_align.axt  -o  japonica_indica_result.txt 
cut -f 1,2,3,4,5 japonica_indica_result.txt  > japonica_indica_result_KaKs.txt

#---------------------------------indica-wheat示例----------------------------------

awk '/^>/ {print $1; next} {print}' wheat.pep.fasta > formatted_wheat.pep.fasta
awk '/^>/ {print $1; next} {print}' wheat.cds.fasta > formatted_wheat.cds.fasta

cat formatted_indica.cds.fasta formatted_wheat.cds.fasta > input_indica_wheat_cds.fasta
cat formatted_indica.pep.fasta formatted_wheat.pep.fasta > input_indica_wheat_pep.fasta

## 数据量比较大，提取前100对ortholog进行分分析
head -n 10 indica_wheat_homo_pairs.txt > indica_wheat_homo_pairs.txt1

## 蛋白比对转cds
echo "6" > proc.txt
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  ParaAT_mdf.pl -h indica_wheat_homo_pairs.txt  -a input_indica_wheat_pep.fasta -n input_indica_wheat_cds.fasta -p proc.txt  -o  indica_wheat_align_out  -m muscle -f axt

cat indica_wheat_align_out/*.axt > indica_wheat_merge_align.axt

## calculate kaks 
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif KaKs_Calculator -m YN -i indica_wheat_merge_align.axt  -o  indica_wheat_result.txt 


#---------------------------------wheat-maize示例----------------------------------
sed -i 's/_P00/_T00/g' maize.pep.fasta
awk '/^>/ {print $1; next} {print}' maize.pep.fasta > formatted_maize.pep.fasta
awk '/^>/ {print $1; next} {print}' maize.cds.fasta > formatted_maize.cds.fasta

cat formatted_wheat.cds.fasta formatted_maize.cds.fasta > input_wheat_maize_cds.fasta
cat formatted_wheat.pep.fasta formatted_maize.pep.fasta > input_wheat_maize_pep.fasta

## 数据量比较大，提取前100对ortholog进行分分析
head -n 10 wheat_maize_homo_pairs.txt > indica_wheat_maize_homo_pairs.txt1

## 蛋白比对转cds
echo "6" > proc.txt
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  ParaAT_mdf.pl -h wheat_maize_homo_pairs.txt  -a input_wheat_maize_pep.fasta -n input_wheat_maize_cds.fasta -p proc.txt  -o  wheat_maize_align_out  -m muscle -f axt

cat wheat_maize_align_out/*.axt > wheat_maize_merge_align.axt

## calculate kaks 
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif KaKs_Calculator -m YN -i wheat_maize_merge_align.axt  -o  wheat_maize_result.txt 

#---------------------------------maize-sorghum示例----------------------------------

awk '/^>/ {print $1; next} {print}' sorghum.pep.fasta > formatted_sorghum.pep.fasta
awk '/^>/ {print $1; next} {print}' sorghum.cds.fasta > formatted_sorghum.cds.fasta

cat formatted_maize.cds.fasta formatted_sorghum.cds.fasta > input_maize_sorghum_cds.fasta
cat formatted_maize.pep.fasta formatted_sorghum.pep.fasta > input_maize_sorghum_pep.fasta

## 数据量比较大，提取前100对ortholog进行分分析
head -n 10 maize_sorghum_homo_pairs.txt > indica_maize_sorghum_homo_pairs.txt1

## 蛋白比对转cds
echo "6" > proc.txt
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif  ParaAT_mdf.pl -h maize_sorghum_homo_pairs.txt  -a input_maize_sorghum_pep.fasta -n input_maize_sorghum_cds.fasta -p proc.txt  -o  maize_sorghum_align_out  -m muscle -f axt

cat maize_sorghum_align_out/*.axt > maize_sorghum_merge_align.axt

## calculate kaks 
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/PhyloTools.sif KaKs_Calculator -m YN -i maize_sorghum_merge_align.axt  -o  maize_sorghum_result.txt 

#==========================================WGD-ortholog 物种间=========================================

#--------------------------------------japonica_indica------------------------------------
cds1=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_japonica.cds.fasta
cds2=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_indica.cds.fasta
outpre=japonica_indica

ln -s  $cds1  input1.cds.fasta
ln -s  $cds2  input2.cds.fasta

## diamond 比对,提取直系同源基因
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  dmd  --nostrictcds -e 1e-10 -o $outpre.dmd  input1.cds.fasta  input2.cds.fasta 

## 计算Ks值
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  ksd --n_threads 10  -o $outpre.ksd  $outpre.dmd/input1.cds.fasta_input2.cds.fasta.rbh  input1.cds.fasta input2.cds.fasta

# head -n 500 $outpre.dmd/input1.cds.fasta_input2.cds.fasta.rbh >  $outpre.dmd/input1.cds.fasta_input2.cds.fasta.rbh500
# singularity  exec  ../software/wgd.sif  wgd  ksd --n_threads 10  -o $outpre.ksd  $outpre.dmd/input1.cds.fasta_input2.cds.fasta.rbh500  input1.cds.fasta input2.cds.fasta

#--------------------------------------indica_wheat------------------------------------
rm input1.cds.fasta
rm input2.cds.fasta

cds1=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_indica.cds.fasta
cds2=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_wheat.cds.fasta
outpre=indica_wheat

ln -s  $cds1  input1.cds.fasta
ln -s  $cds2  input2.cds.fasta

## diamond 比对,提取直系同源基因
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  dmd  --nostrictcds -e 1e-10 -o $outpre.dmd  input1.cds.fasta  input2.cds.fasta 

## 计算Ks值
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  ksd --n_threads 10  -o $outpre.ksd  $outpre.dmd/input1.cds.fasta_input2.cds.fasta.rbh  input1.cds.fasta input2.cds.fasta


#--------------------------------------wheat_maize------------------------------------
rm input1.cds.fasta
rm input2.cds.fasta
cds1=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_wheat.cds.fasta
cds2=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_maize.cds.fasta
outpre=wheat_maize

ln -s  $cds1  input1.cds.fasta
ln -s  $cds2  input2.cds.fasta

## diamond 比对,提取直系同源基因
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  dmd  --nostrictcds -e 1e-10 -o $outpre.dmd  input1.cds.fasta  input2.cds.fasta 

## 计算Ks值
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  ksd --n_threads 10  -o $outpre.ksd  $outpre.dmd/input1.cds.fasta_input2.cds.fasta.rbh  input1.cds.fasta input2.cds.fasta


#--------------------------------------maize_sorghum------------------------------------
rm input1.cds.fasta
rm input2.cds.fasta
cds1=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_maize.cds.fasta
cds2=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_sorghum.cds.fasta
outpre=maize_sorghum

ln -s  $cds1  input1.cds.fasta
ln -s  $cds2  input2.cds.fasta

## diamond 比对,提取直系同源基因
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  dmd  --nostrictcds -e 1e-10 -o $outpre.dmd  input1.cds.fasta  input2.cds.fasta 

## 计算Ks值
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  ksd --n_threads 10  -o $outpre.ksd  $outpre.dmd/input1.cds.fasta_input2.cds.fasta.rbh  input1.cds.fasta input2.cds.fasta


#==========================================WGD-paralog 物种内部=========================================
#-------------------------------------indica--------------------------------------
# ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/6Collinearity/MCScanx/indica.gff ./indica.gff
# ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/6Collinearity/MCScanx/japonica.gff ./japonica.gff
# ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/6Collinearity/MCScanx/wheat.gff ./wheat.gff
# ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/6Collinearity/MCScanx/maize.gff ./maize.gff
# ln -s /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/6Collinearity/MCScanx/sorghum.gff ./sorghum.gff



cds=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_indica.cds.fasta
outpre=indica

ln -s  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_indica.cds.fasta  input.cds.fasta

## diamond 比对 并进行MCL聚类
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd dmd  -e 1e-10  --nostrictcds -o $outpre.dmd  input.cds.fasta 

## 计算Ks值
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  ksd --n_threads 10 -mp 1000 -o $outpre.ksd  $outpre.dmd/input.cds.fasta.mcl input.cds.fasta

## 共线性分析,提取共线性基因对Ks值 
#singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd syn -o $outpre.syn  -ks $outpre.ksd/input.cds.fasta.ks.tsv  -f mRNA -a ID input.gff  $outpre.dmd/input.cds.fasta.mcl

# ## 基于mcscanx结果,提取共线性基因对Ks值 
cp  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/6Collinearity/MCScanx/indica.collinearity  ./
perl /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/filter_ks_mcscanx.pl   indica.ksd/input.cds.fasta.ks.tsv  indica.collinearity > indica.syn.ks.tsv
#-------------------------------------japonica--------------------------------------
rm input.cds.fasta
cds=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_japonica.cds.fasta
outpre=japonica

ln -s  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_japonica.cds.fasta  input.cds.fasta

## diamond 比对 并进行MCL聚类
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd dmd  -e 1e-10  --nostrictcds -o $outpre.dmd  input.cds.fasta 

## 计算Ks值
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  ksd --n_threads 10 -mp 1000 -o $outpre.ksd  $outpre.dmd/input.cds.fasta.mcl input.cds.fasta

## 共线性分析,提取共线性基因对Ks值 
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd syn -o $outpre.syn  -ks $outpre.ksd/input.cds.fasta.ks.tsv  -f mRNA -a ID input.gff  $outpre.dmd/input.cds.fasta.mcl

#-------------------------------------wheat--------------------------------------
rm input.cds.fasta
cds=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_wheat.cds.fasta
outpre=wheat

ln -s  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_wheat.cds.fasta  input.cds.fasta

## diamond 比对 并进行MCL聚类
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd dmd  -e 1e-10  --nostrictcds -o $outpre.dmd  input.cds.fasta 

## 计算Ks值
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  ksd --n_threads 10 -mp 1000 -o $outpre.ksd  $outpre.dmd/input.cds.fasta.mcl input.cds.fasta

## 共线性分析,提取共线性基因对Ks值 
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd syn -o $outpre.syn  -ks $outpre.ksd/input.cds.fasta.ks.tsv  -f mRNA -a ID input.gff  $outpre.dmd/input.cds.fasta.mcl

#-------------------------------------maize--------------------------------------
rm input.cds.fasta
cds=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_maize.cds.fasta
outpre=maize

ln -s  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_maize.cds.fasta  input.cds.fasta

## diamond 比对 并进行MCL聚类
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd dmd  -e 1e-10  --nostrictcds -o $outpre.dmd  input.cds.fasta 

## 计算Ks值
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  ksd --n_threads 10 -mp 1000 -o $outpre.ksd  $outpre.dmd/input.cds.fasta.mcl input.cds.fasta

## 共线性分析,提取共线性基因对Ks值 
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd syn -o $outpre.syn  -ks $outpre.ksd/input.cds.fasta.ks.tsv  -f mRNA -a ID input.gff  $outpre.dmd/input.cds.fasta.mcl

#-------------------------------------sorghum--------------------------------------

rm input.cds.fasta
cds=/home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_sorghum.cds.fasta
outpre=sorghum

ln -s  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/7Kaks/formatted_sorghum.cds.fasta  input.cds.fasta

## diamond 比对 并进行MCL聚类
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd dmd  -e 1e-10  --nostrictcds -o $outpre.dmd  input.cds.fasta 

## 计算Ks值
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd  ksd --n_threads 10 -mp 1000 -o $outpre.ksd  $outpre.dmd/input.cds.fasta.mcl input.cds.fasta

## 共线性分析,提取共线性基因对Ks值 
singularity  exec  /home/win/4T/GeneDL/DXZ_DL/expriments/CompareGenome/Comparative_genomics/software/wgd.sif  wgd syn -o $outpre.syn  -ks $outpre.ksd/input.cds.fasta.ks.tsv  -f mRNA -a ID input.gff  $outpre.dmd/input.cds.fasta.mcl

## 基于mcscanx结果,提取共线性基因对Ks值 
cp  ../05.Collinearity/MCScanX/sind.collinearity  ./
perl ../software/filter_ks_mcscanx.pl   Sind.ksd/input.cds.fasta.ks.tsv  sind.collinearity > Sind.syn.ks.tsv

#=================================可以采用单物种，也可以用多物种======================================
#circos绘图


#LTR分析


