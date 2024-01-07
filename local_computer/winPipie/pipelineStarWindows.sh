# STAR转录组比对和定量 {#STAR}

## 服务器登录

## IP: 192.168.1.107 端口 22

## 用户名: 姓名汉语拼音

## 密码: yishengxin

# 设置工作目录
wd=/mnt/d/train/serverData/data
  
# 激活环境
conda activate transcriptome

cd ${wd}

ls -sh *.fq.gz


### 4.0K compare_pair          25M trt_N080611_1.fq.gz     16M untrt_N061011_1.fq.gz
### 4.0K sampleFile            25M trt_N080611_2.fq.gz     16M untrt_N061011_2.fq.gz
###  13M trt_N052611_1.fq.gz   14M trt_N61311_1.fq.gz      19M untrt_N080611_1.fq.gz
###  13M trt_N052611_2.fq.gz   14M trt_N61311_2.fq.gz      19M untrt_N080611_2.fq.gz
###  18M trt_N061011_1.fq.gz   19M untrt_N052611_1.fq.gz   15M untrt_N61311_1.fq.gz
###  18M trt_N061011_2.fq.gz   19M untrt_N052611_2.fq.gz   15M untrt_N61311_2.fq.gz
### 
### genome:
### 总用量 87M
###  62M GRCh38.fa   24M GRCh38.gtf  1.7M GRCh38.idmap


## 测序文件下载

# 使用NCBI提供的SRA-toolkit中的工具`fastq-dump`直接下载SRR文件，并转换为`FASTQ`格式，`--split-3`参数表示如果是双端测序就自动拆分，如果是单端不受影响。

# SRA toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>, 根据服务器操作系统类型下载对应的二进制编码包，下载解压放到环境变量即可使用。

# CentOS下地址：<https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.0/sratoolkit.2.9.0-centos_linux64.tar.gz>。

# 如果需要下载测序数据，则去掉下面语句前的#
#cd ${wd}
#fastq-dump -v --split-3 --gzip SRR1039508
#rename "SRR1039508"  "untrt_N61311"  SRR1039508*
#fastq-dump -v --split-3 --gzip SRR1039509
#rename "SRR1039509"  "trt_N61311"  SRR1039509*
#fastq-dump -v --split-3 --gzip SRR1039512
#rename "SRR1039512"  "untrt_N052611"  SRR1039512*
#fastq-dump -v --split-3 --gzip SRR1039513
#rename "SRR1039513"  "trt_N052611"  SRR1039513*
#fastq-dump -v --split-3 --gzip SRR1039516
#rename "SRR1039516"  "untrt_N080611"  SRR1039516*
#fastq-dump -v --split-3 --gzip SRR1039517
#rename "SRR1039517"  "trt_N080611"  SRR1039517*
#fastq-dump -v --split-3 --gzip SRR1039520
#rename "SRR1039520"  "untrt_N061011"  SRR1039520*
#fastq-dump -v --split-3 --gzip SRR1039521
#rename "SRR1039521"  "trt_N061011"  SRR1039521*
#rename "fastq" "fq" *.gz
#/bin/rm ~/ncbi/public/sra/*.sra

## 测序文件查看和测序量评估

# 进入data目录
cd ${wd}

# 查看文件，压缩或未压缩的都可以，尤其适合查看极大文件
# 在Rstudio界面不可用, 可以查看，但不能退出
#less trt_N061011_1.fq.gz

# zcat查看gzip压缩的文件
# head -n 8 显示前8行文件内容
# 前8行是代表几条序列？
zcat trt_N061011_1.fq.gz | head -n 8

### @SRR1039521.13952745/1
### TTCCTTCCTCCTCTCCCTCCCTCCCTCCTTTCTTTCTTCCTGTGGTTTTTTCCTCTCTTCTTC
### +
### HIJIIJHGHHIJIIIJJJJJJJJJJJJJJJJJJJJJIIJJFIDHIBGHJIHHHHHHFFFFFFE
### @SRR1039521.7213571/1
### GAGGAAGGGCAGAGGGAGCAGGGAGACTGTAGATCAGGGGCTGAATGGAGATCCGGTCCTG
### +
### <1E:E@3B:8?E?B;6FFFI@@6/0@-B<FFD4)=C=D>'-54@B@>>C>AC;>=38=6A@

# 测序reads数计算

# wc -l: 计算行数
# bc -l: 计算器 (-l：浮点运算)
# 为什么除以4，又除以1000000？
echo "`zcat trt_N061011_1.fq.gz | wc -l` / (4*1000000)" | bc -l
# .464329 million

# 测序碱基数计算 (只包含20号染色体的数据) (awk的介绍见[常用和不太常用的awk命令](http://mp.weixin.qq.com/s/8wD14FXt7fLDo1BjJyT0ew))

# awk运算
# %取余数
zcat trt_N061011_1.fq.gz | awk '{if(FNR%4==0) base+=length}END{print base/10^9,"G";}'
# 0.028816 G

## 测序质量评估 (NGS基础 - FASTQ格式解释和质量评估: https://mp.weixin.qq.com/s/tDMih7ISLJcL4F4sWBq3Vw)

# fastqc trt_N061011_1.fq.gz
# trt_N061011_1_fastqc.html 在Rstudio中打开 (View in Web browser)

# 批量评估
# fastqc *.fq.gz

# multiqc 整理评估结果
# multiqc -d . -o multiqc

# 结果存储在multiqc/multiqc_report.html, 双击打开

## 基因注释文件准备 gtf转bed12
cd ${wd}/genome

# gtf2bed12.sh所做的操作就是下面两句话
# gtfToGenePred -ignoreGroupsWithoutExons GRCh38.gtf GRCh38.gtf.50505050.pred
# genePredToBed GRCh38.gtf.50505050.pred GRCh38.gtf.bed12
# 删除临时文件
# /bin/rm -f GRCh38.gtf.50505050.pred

# 
gtf2bed12.sh -f GRCh38.gtf
# 选择长度适中的转录本用于后续评估
awk '$3-$2>1000 && $3-$2<2000' GRCh38.gtf.bed12 >GRCh38.model.gtf.bed12
# head GRCh38.gtf.bed12

# 生成igv截图脚本
cut -f 1-6 GRCh38.gtf.bed12 >igv.capture.bed
bedtools igv -i igv.capture.bed -path "C:" >igv.batch.script

# # cp ../bak/GRCh38.chromsize .
# 下载chrome size文件，第一列为染色体名字，第二列为染色体大小
# 也可以自己写程序计算
# 适用于任何基因组序列，先判断是不是读到的第一条染色体，如果不是就打印；如果是只存储
# 记得输出最后一条染色体的长度信息
# awk 'BEGIN{OFS="\t"}{if($0~/>/) {if(size>0) print chrname, size; size=0; chrname=$0; sub(">","", chrname);} else size+=length;}END{print chrname,size}' GRCh38.fa
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg38.chromInfo" | tail -n +2 >GRCh38.chromsize
# head GRCh38.chromsize


# STAR构建基因组索引 (基因组索引构建的目的是为了方便快速比对)
# 需要GTF文件中有第三列的区域类型中有exon，第9列的区域属性信息中有gene_id\transcript_id
cd ${wd}/genome

mkdir -p star_GRCh38
# --runThreadN 2: 指定使用2个线程
# --sjdbOverhang 100: 默认

# 如果代码被分割为多行（如下，行末有 \ ），建议全选运行

STAR --runMode genomeGenerate --runThreadN 2 --genomeDir star_GRCh38 \
     --genomeFastaFiles GRCh38.fa --sjdbGTFfile GRCh38.gtf &


# Aug 03 10:52:41 ..... started STAR run
# Aug 03 10:52:41 ... starting to generate Genome files
# Aug 03 10:52:42 ... starting to sort Suffix Array. This may take a long time...
# Aug 03 10:52:43 ... sorting Suffix Array chunks and saving them to disk...
# Aug 03 10:53:01 ... loading chunks from disk, packing SA...
# Aug 03 10:53:04 ... finished generating suffix array
# Aug 03 10:53:04 ... generating Suffix Array index
# Aug 03 10:53:33 ... completed Suffix Array index
# Aug 03 10:53:33 ..... processing annotations GTF
# Aug 03 10:53:33 ..... inserting junctions into the genome indices
# Aug 03 10:53:50 ... writing Genome to disk ...
# Aug 03 10:53:50 ... writing Suffix Array to disk ...
# Aug 03 10:53:50 ... writing SAindex to disk
# Aug 03 10:53:51 ..... finished successfully

# 输出文件如下，注意看下输出文件的大小，有无空文件，基因数是否对。

ls -sh star_GRCh38

# 总用量 2.1G
# 4.0K chrLength.txt      368K exonInfo.tab          1.5G SAindex
# 4.0K chrNameLength.txt   24K geneInfo.tab          204K sjdbInfo.txt
# 4.0K chrName.txt         64M Genome                204K sjdbList.fromGTF.out.tab
# 4.0K chrStart.txt       4.0K genomeParameters.txt  204K sjdbList.out.tab
# 732K exonGeTrInfo.tab   516M SA                    224K transcriptInfo.tab

cp star_GRCh38/chrNameLength.txt GRCh38.chromsize

head GRCh38.chromsize

# STAR解析后的基因数
head -n 1 star_GRCh38/geneInfo.tab

# 原始GTF的基因数
grep -cP '\tgene\t' GRCh38.gtf

## Reads比对

# mkdir新建目录
cd ${wd}
mkdir -p trt_N061011
# --runThreadN 4: 使用4个线程
# --readFilesIn: 输入文件，左端和右端
# --readFilesCommand zcat：gzip压缩，指定解压方式
# --genomeDir：基因组索引目录的位置
# -S: 指定输出文件

# 动物一般写 1000000，植物一般写5000
max_intron_size=1000000

# --genomeLoad LoadAndKeep : 共享内存
# 其他参数自己对着star的帮助手册看

# ***这里不要忘记运行**

# star_p=" --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
#       --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
#       --alignIntronMin 20 --alignIntronMax ${max_intron_size} \
#                --alignMatesGapMax ${max_intron_size} \
#                --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
#                --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
#                --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
#                --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
#                --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts"

# STAR比对单个样品 
# # 动物一般写 1000000，植物一般写5000
# 对应参数 alignIntronMax 和 alignMatesGapMax
# max_intron_size=1000000
STAR --runMode alignReads --runThreadN 4 \
        --readFilesIn trt_N061011_1.fq.gz trt_N061011_2.fq.gz \
        --readFilesCommand zcat --genomeDir genome/star_GRCh38 \
        --outFileNamePrefix trt_N061011/trt_N061011.  --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
       --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
       --alignIntronMin 20 --alignIntronMax 1000000 \
       --alignMatesGapMax 1000000 \
       --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
       --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
       --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
       --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
       --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts \
       --outTmpDir /tmp/trt_N061011/


# Aug 03 11:44:27 ..... started STAR run
# Aug 03 11:44:27 ..... loading genome
# Aug 03 11:44:30 ..... started mapping
# Aug 03 11:44:48 ..... finished successfully

# 输出文件

ls -sh trt_N061011/*

# trt_N061011.Aligned.out.bam: 比对到基因组的bam文件
# trt_N061011.Aligned.toTranscriptome.out.bam：比对到转录组的bam文件，供RSEM计算TPM使用 
# trt_N061011.Log.final.out: reads比对日志
# trt_N061011.Log.out: 运行参数和过程
# trt_N061011.Log.progress.out: 运行日志

# trt_N061011.ReadsPerGene.out.tab: 每个基因的reads count，链非特异性RNASeq选第2列.
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
# 
# Select the output according to the strandedness of your data. Note, that if you have stranded data and choose one of the columns 3 or 4, the other column (4 or 3) will give you the count of antisense reads. 

# trt_N061011.SJ.out.tab: Junction reads
# trt_N061011.Unmapped.out.mate1：未比对上的reads
# trt_N061011.Unmapped.out.mate2：未比对上的reads

# using the --outFilterType Normal vs BySJout options will create slightly different SJ.out.tab counts for the following reason.
# Imagine that you have a read two junctions, with only 1st junction passing the filter.
# Then with the Normal option, the 1st junction will be counted in the SJ.out.tab.
# If the BySJout option is used, the entire read alignment may be prohibited because of the 2nd junction, and then the 1st junction will be not counted in the SJ.out.tab.

#筛选reads，按坐标排序、索引BAM文件供下游使用, 也可导入IGV查看reads比对情况、是否有变异位点等。

# samtools具体参数解释见 samtools -?
# -@ 4: 4个线程
cd ${wd}
mkdir -p tmp
samtools sort -@ 4 -T tmp/trt_N061011 \
  -o trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
  trt_N061011/trt_N061011.Aligned.out.bam
samtools index trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam
# 为什么要按坐标排序？
# 为什么要建索引？
# 就可以导入IGV中查看reads的比对情况了

# BigWig峰图文件生成，导入IGV或UCSC genomebrowser获取表达丰度图。

# Wig里面有什么？
# 为什么要生成BigWig？
# 是否需要标准化？

cd ${wd}
STAR --runMode inputAlignmentsFromBAM \
        --inputBAMfile trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
        --outWigType bedGraph --outFileNamePrefix trt_N061011/trt_N061011. \
        --outWigNorm RPM --outWigStrand Unstranded
bedSort trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bg \
        trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bg
bedGraphToBigWig trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bg \
        genome/star_GRCh38/chrNameLength.txt \
        trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bw

# bam2wig.py是RseQC中的软件，可以使用`which bam2wig.py`查看命令的路径
# -s: 是chrom size文件，两列文件，前面又介绍
# -t: 归一化因子，所有样品都归一化到统一测序深度，方便比较
# bam2wig.py -i trt_N061011/trt_N061011.sortP.bam -s genome/GRCh38.chromsize -o trt_N061011/trt_N061011 -t 1000000000 -q 0
# ls -ltr trt_N061011/trt_N061011*

## MultiQC查看STAR比对结果

multiqc -f -d . -o multiqc

### 比对质量评估
#
#### Reads在基因上的分布评估
#
## 程序运行结束后，默认生成折线图
## 如果样品多，也可以生成理论课件中的热图
# 如果缺少libreadline.so.6则执行 ln -s ~/anaconda3/envs/transcriptome/lib64/libreadline.so ~/anaconda3/envs/transcriptome/lib64/libreadline.so.6
# 如果缺少libncurses.so.5则执行 ln -s ~/anaconda3/envs/transcriptome/lib64/libncurses.so ~/anaconda3/envs/transcriptome/lib64/libncurses.so.5
geneBody_coverage2.py -i \
 trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bw \
 -r genome/GRCh38.model.gtf.bed12 -o trt_N061011/trt_N061011.geneBody_coverage
#
#### Reads比对到基因组标志性区域的分布
#
## 请完成堆积柱状图的绘制
## 或使用www.ehbio.com/ImageGP
read_distribution.py -i trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
  -r genome/GRCh38.gtf.bed12 >trt_N061011/trt_N061011.read_distrib.xls
#
cat trt_N061011/trt_N061011.read_distrib.xls
#
#
#### 测序饱和度评估
#
## -s: 采样频率，0-100之间的整数，类似于步长
## -q: 过滤低质量比对
RPKM_saturation.py -i \
 trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
 -r genome/GRCh38.gtf.bed12 -s 10 -q 0 -o trt_N061011/trt_N061011.RPKM_saturation
#
#ls -ltr trt_N061011
#---------------------------------------------------------

## 批量比对和定量

##先自己建立样品分组信息文件，命名为`sampleFile` (可以是任意名字, 后面做差异分析也会用到)

## 内容如下：第一行为标题行；第一列为样品名字；后面为样品信息；此处只用到第一列。第二列名字尽量为`conditions`，方便后期直接进行差异基因检测（如果您懂得原理，可以写任何列名字）。

#   Samp    conditions
#   untrt_N61311    untrt
#   untrt_N052611   untrt
#   untrt_N080611   untrt
#   untrt_N061011   untrt
#   trt_N61311      trt
#   trt_N052611     trt
#   trt_N080611     trt
#   trt_N061011     trt

# TAB键分割
#cat <<END | sed 's/|/\t/' >sampleFile
#Samp|conditions
#untrt_N61311|untrt
#untrt_N052611|untrt
#untrt_N080611|untrt
#untrt_N061011|untrt
#trt_N61311|trt
#trt_N052611|trt
#trt_N080611|trt
#trt_N061011|trt
#END


# 测试是否为TAB键分割
cat -A sampleFile
# Samp^Iconditions$
# untrt_N61311^Iuntrt$
# untrt_N052611^Iuntrt$
# untrt_N080611^Iuntrt$
# untrt_N061011^Iuntrt$
# trt_N61311^Itrt$
# trt_N052611^Itrt$
# trt_N080611^Itrt$
# trt_N061011^Itrt$

cd ${wd}

# 动物一般写 1000000，植物一般写5000
# max_intron_size=1000000

# --genomeLoad LoadAndKeep : 共享内存
# 其他参数自己对着star的帮助手册看
# star_p=" --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
#       --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
#       --alignIntronMin 20 --alignIntronMax ${max_intron_size} \
#                --alignMatesGapMax ${max_intron_size} \
#                --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
#                --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
#                --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
#                --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
#                --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts"

for i in `tail -n +2 sampleFile | cut -f 1`; do 
	mkdir -p ${i}
	mkdir -p tmp
	STAR --runMode alignReads --runThreadN 4 \
        --readFilesIn ${i}_1.fq.gz ${i}_2.fq.gz \
        --readFilesCommand zcat --genomeDir genome/star_GRCh38 \
        --outFileNamePrefix ${i}/${i}. --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
       --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
       --alignIntronMin 20 --alignIntronMax 1000000 \
       --alignMatesGapMax 1000000 \
       --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
       --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
       --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
       --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
       --outTmpDir /tmp/${i}/ \
       --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts
  samtools sort -@ 10 -T ${i}.tmp \
                -o ${i}/${i}.Aligned.sortedByCoord.out.bam \
                ${i}/${i}.Aligned.out.bam
	samtools index ${i}/${i}.Aligned.sortedByCoord.out.bam
	STAR --runMode inputAlignmentsFromBAM \
        --inputBAMfile ${i}/${i}.Aligned.sortedByCoord.out.bam \
        --outWigType bedGraph --outFileNamePrefix ${i}/${i}. \
        --outWigNorm RPM --outWigStrand Unstranded
  bedSort ${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
        ${i}/${i}.Signal.UniqueMultiple.str1.out.bg
  bedGraphToBigWig ${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
        genome/star_GRCh38/chrNameLength.txt \
        ${i}/${i}.Signal.UniqueMultiple.str1.out.bw
done &



## 批量评估

# cd ${wd}
# for i in `tail -n +2 sampleFile | cut -f 1`; do 
#   /anaconda2/bin/geneBody_coverage2.py -i \
#     ${i}/${i}.Signal.UniqueMultiple.str1.out.bw \
#     -r genome/GRCh38.model.gtf.bed12 -o ${i}/${i}.geneBody_coverage
#   /anaconda2/bin/read_distribution.py -i ${i}/${i}.Aligned.sortedByCoord.out.bam \
#     -r genome/GRCh38.gtf.bed12 >${i}/${i}.read_distrib.xls
#   /anaconda2/bin/RPKM_saturation.py -i \
#     ${i}/${i}.Aligned.sortedByCoord.out.bam \
#     -r genome/GRCh38.gtf.bed12 -s 10 -q 0 -o ${i}/${i}.RPKM_saturation
# done

## multiqc整理软件运行结果

multiqc -f -d . -o multiqc

## 合并表达文件

# 基因reads count增加样品信息，方便后续合并        
# sed '5 i\Gene\ttrt_N061011\ttrt_N061011\ttrt_N061011' trt_N061011/trt_N061011.ReadsPerGene.out.tab trt_N061011/trt_N061011.ReadsPerGene.out.tab.ehbio

cd ${wd}
for i in `tail -n +2 sampleFile | cut -f 1`; do 
  sed "5 i\Gene\t${i}\t${i}\t${i}" ${i}/${i}.ReadsPerGene.out.tab >${i}/${i}.ReadsPerGene.out.tab.ehbio
done


# Linux 命令合并
cd ${wd}
paste `find . -name *.ReadsPerGene.out.tab.ehbio` | tail -n +5 | \
  awk 'BEGIN{OFS=FS="\t" }{line=$1; \
    for(i=2;i<=NF;i++) if(i%2==0 && i%4!=0) line=line"\t"$i; print line;}' \
  >ehbio_trans.Count_matrix.xls
head ehbio_trans.Count_matrix.xls



## 番外：multiqc重命名

# General statics 重命名
for i in `tail -n +2 sampleFile | cut -f 1`; do echo -e "$i | $i\tSTAR_${i}"; done >star_map
for i in `tail -n +2 sampleFile | cut -f 1`; do echo -e "$i | $i.salmon.count | aux_info | ${i}.salmon.count\tSalmon_${i}"; done >salmon_map
# 在Rstudio中打开这两个文件，拷贝到MultiQC界面
#




