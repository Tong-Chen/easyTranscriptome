# Salmon定量 {#salmon}

	# 1. 这一部分是在远程服务器的Rstudio界面利用远程服务器的计算资源进行的操作。
	# 2. 在服务器运行时，不能跳着运行，有时前一步的输出结果是后一步的输入结果。
	# 3. 文件的上传和下载演示。
	# 4. 此服务器只供练习使用，有效期一个月。
	# 5. 所有#号开头的行是注释
	# 6. 此文为Markdown格式，放入Markdown阅读器可查看其层级结构

## 服务器登录

#	IP: 192.168.1.107 端口 22

#	用户名: 姓名汉语拼音

#	密码: yishengxin

  wd=/mnt/d/train/serverData/data
  
  # 激活环境
  conda activate transcriptome

## 切换目录，查看文件

	cd ${wd}

	# -sh 表示获得每个文件的大小(s: size)，以人可读的形式显示 (h)即加上单位如M，K，G等。
	ls -sh *

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


## 测序原始数据下载

	# 使用NCBI提供的SRA-toolkit中的工具`fastq-dump`直接下载SRR文件，
	# 并转换为`FASTQ`格式，`--split-3`参数表示如果是双端测序就自动拆分，如果是单端不受影响。

	# SRA toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>, 
	# 根据服务器操作系统类型下载对应的二进制编码包，下载解压放到环境变量即可使用。

	# CentOS下地址：<https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.0/sratoolkit.2.9.0-centos_linux64.tar.gz>。

	# 如果需要下载测序数据，则取消下面语句前的#
	# 参加线上课程的老师，请取消下面的注释，自行下载

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

## 测序原始数据查看和测序量评估

	# 输入文件：测序公司返回的FASTQ文件 *.fq.gz 如trt_N061011_1.fq.gz
	# 输出文件：html文件 如 trt_N061011_1_fastqc.html

	## 测序质量评估 (NGS基础 - FASTQ格式解释和质量评估: https://mp.weixin.qq.com/s/tDMih7ISLJcL4F4sWBq3Vw)
	cd ${wd}
	fastqc trt_N061011_1.fq.gz
	# trt_N061011_1_fastqc.html 在Rstudio中打开 (View in Web browser)

	# 批量评估
	fastqc *.fq.gz

	# multiqc 整理评估结果
	# -d .: 表示分析当前目录所有文件和文件夹，multiqc会遍历读取，分析每一个文件是不是
	#       常用工具的输出结果文件，再进行读取
	multiqc -d . -o multiqc

	# 结果存储在multiqc/multiqc_report.html, 双击打开


# 不基于比对的定量 (salmon)

## 基因注释文件准备

	# 工作目录：${wd}/genome
	# 输入文件：
	#		GRCh38.gtf 人基因注释序列，从Ensembl下载
	#       GRCh38.fa 人基因组序列，从Ensembl下载
	# 输出文件：
	#       GRCh38.transcript.fa 所有转录本序列文件

	cd ${wd}/genome

	## 获取转录本序列
	
	# 如果物种有cDNA序列提供，直接下载；但Biomart直接下载的cDNA只包含蛋白编码基因，建议自己提取
	
	# 如果是三代测序，可直接用公司提供的序列作为cDNA
	
	# 如果是无参转录组，直接用Trinity拼装的转录本作为cDNA
	
	# 如果没有则用下面的程序，通过gtf和genome文件提取

	# GRCh38.fa 人基因组序列，从Ensembl下载
	# GRCh38.gtf 人基因注释序列，从Ensembl下载
	# 
	# wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz -O GRCh38.fa.gz
	# gunzip -c GRCh38.fa.gz >GRCh38.fa
	# wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O GRCh38.gtf.gz
	# gunzip -c GRCh38.gtf.gz >GRCh38.gtf

	gffread GRCh38.gtf -g GRCh38.fa -w GRCh38.transcript.fa.tmp

	# gffread生成的fasta文件同时包含基因名字和转录本名字
	grep '>' GRCh38.transcript.fa.tmp | head

	# 去掉空格后面的字符串，保证cDNA文件中fasta序列的名字简洁，不然后续会出错
	cut -f 1 -d ' ' GRCh38.transcript.fa.tmp >GRCh38.transcript.fa
	head GRCh38.transcript.fa

	grep '>' GRCh38.transcript.fa | head

## 构建Salmon索引 (cDNA)

	# 工作目录：${wd}/genome
	# 输入文件：
	#       GRCh38.transcript.fa 所有转录本序列文件
	# 输出文件：
	#		GRCh38.salmon：一个文件夹，包含salmon的索引
			
	# GRCh38.transcript.fa：待分析的cDNA文件，转录本序列fasta文件
	# GRCh38.salmon： salmon索引输出文件夹
	salmon index -t GRCh38.transcript.fa -i GRCh38.salmon

	# salmon索引目录下的内容
	ls GRCh38.salmon

  # complete_ref_lens.bin  duplicate_clusters.tsv  pos.bin           refAccumLengths.bin  refseq.bin
  # ctable.bin             info.json               pre_indexing.log  ref_indexing.log     seq.bin
  # ctg_offsets.bin        mphf.bin                rank.bin          reflengths.bin       versionInfo.json

	# 有些转录本名字不同，但序列完全一样，被鉴定出来只保留一个。
	head GRCh38.salmon/duplicate_clusters.tsv 

	# RetainedTxp	DuplicateTxp
	# ENST00000486475	ENST00000629881
	# ENST00000401322	ENST00000637366
	# ENST00000401322	ENST00000401376
	# ENST00000401322	ENST00000621667

## 构建Salmon索引 (cDNA+genome)

  # 区分序列是来源于转录本还是其它基因组区域

  # 参考
  # https://github.com/COMBINE-lab/salmon
  # https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

  # 获取所有基因组序列的名字存储于decoy中
  grep '^>' GRCh38.fa | cut -d ' ' -f 1 | sed 's/^>//g' >GRCh38.decoys.txt

  # 合并cDNA和基因组序列一起
  # 注意：cDNA在前，基因组在后
  
  cat GRCh38.transcript.fa GRCh38.fa >GRCh38_trans_genome.fa
  
  # 构建索引 （更慢，结果会更准） # --no-version-check
  salmon index -t GRCh38_trans_genome.fa -d GRCh38.decoys.txt -i GRCh38.salmon_sa_index -p 5

## 单个样品Salmon定量获得基因表达表 

	# 工作目录：${wd}/
	# 输入文件：
	#       trt_N061011_1.fq.gz, trt_N061011_2.fq.gz: 双端测序数据
	#       genome/GRCh38.salmon：基因组索引 (上一步构建获得的)
	# 输出文件：
	#		trt_N061011/trt_N061011.salmon.count: salmon定量结果输出目录
	#       trt_N061011/trt_N061011.salmon.count/quant.sf: 样品转录组定量结果，差异分析时会用到

	# 注意切换目录
	cd ${wd}/

	# -p: 表示若待创建的文件夹已存在则跳过；若不存在，则创建；也可用于创建多层文件夹
	# man mkdir 可查看详细帮助
	mkdir -p trt_N061011

	# -l: 自动判断文库类型，尤其适用于链特异性文库
	# The library type -l should be specified on the command line 
	# before the read files (i.e. the parameters to -1 and -2, or -r). 
	# This is because the contents of the library type flag is used to determine how the reads should be interpreted.

	# --gcBias: 校正测序片段GC含量，获得更准确的转录本定量结果
	# One can simply run Salmon with --gcBias in any case, 
	# as it does not impair quantification for samples without GC bias, 
	# it just takes a few more minutes per sample. 
	# For samples with moderate to high GC bias, correction for this bias at the 
	# fragment level has been shown to reduce isoform quantification errors
	salmon quant --gcBias -l A -1 trt_N061011_1.fq.gz -2 trt_N061011_2.fq.gz  -i genome/GRCh38.salmon_sa_index -o trt_N061011/trt_N061011.salmon.count -p 5

	# 如果是单端数据，只需提供 -r
	# salmon quant --gcBias -l A -r trt_N061011.fq.gz -i genome/GRCh38.salmon -o trt_N061011/trt_N061011.salmon.count -p 4

	# 输出结果存储在 trt_N061011/trt_N061011.salmon.count目录中
	# quant.sf 为转录本表达定量结果，第4列为TPM结果，第5列为reads count
	# quant.genes.sf 为基因表达定量结果
	head -n 30 trt_N061011/trt_N061011.salmon.count/quant.sf | tail

## 所有样本批量定量

	# 工作目录：${wd}/
	# 输入文件：
	#		sampleFile: TAB键分割的样品分组信息，至少2列，样本名字和样本所属分组信息
	#                   注意：样本名字与测序数据文件名前缀一致 **
	#       genome/GRCh38.salmon：基因组索引 (上一步构建获得的)
	# 隐式输入文件：
	#       每个样本质控后的测序数据 (对Salmon来讲，提供raw data和clean data差别不大，因为是基于K-mer的计算)
	# 		sampleFile中样本名是 trt_N061011，加上后缀_1.fq.gz可获得左端测序数据名字，加上后缀_2.fq.gz可获得右端测序数据名字。
	#       trt_N061011_1.fq.gz, trt_N061011_2.fq.gz: 双端测序数据
	# 输出文件：
	#       每个样本一个文件夹
	#		样本名/样本名.salmon.count: salmon定量结果输出目录
	#       样本名/样本名.salmon.count/quant.sf: 样品转录组定量结果，差异分析时会用到
	#		trt_N061011/trt_N061011.salmon.count: salmon定量结果输出目录
	#       trt_N061011/trt_N061011.salmon.count/quant.sf: 样品转录组定量结果，差异分析时会用到

	cd ${wd}

	# TAB键分割的样品分组信息
	# 可以用任何文本编辑工具或者excel导出一个TAB键分割的文件，
	# 第一列是样本名字
	# 后面的列是样本分组或其它熟悉信息，记录越全越好
	head sampleFile

	# 批量定量
	# -g genome/GRCh38.gtf: 同时输出基因的定量结果，默认只有转录本定量结果
	cd ${wd}
	
	# 注意末尾的 &，把程序放入后台，避免被取词工具的复制键 (ctrl+c)干扰。ctrl+c在命令行下是打断程序运行
	
	for samp in `tail -n +2 sampleFile | cut -f 1`; do salmon quant --gcBias -l A -1 ${samp}_1.fq.gz -2 ${samp}_2.fq.gz  -i genome/GRCh38.salmon_sa_index -o ${samp}/${samp}.salmon.count -p 4 >${samp}.salmon.log 2>&1; done &

	# 命令拆解，便于理解
	# tail -n +2 sampleFile
	# tail -n +2 sampleFile | cut -f 1
	# for samp in `tail -n +2 sampleFile | cut -f 1`; do echo ${samp}_1.fq.gz; done


## 利用multiqc整理salmon的输出日志

	# 获得reads比对率信息
	# -f覆盖之前的日志
	multiqc -f -d . -o multiqc/

## 整理Salmon输出用于后续差异基因分析

	# 工作目录：${wd}/
	# 用于差异基因分析的3个文件：
	#		quant.sf.zip: 所有样本的quant.sf文件的集合
	#		salmon.output: 一个文本文件，用于记录每个样本和对应的quant.sf所在的路径，方便做差异基因分析时同时读取quant.sf.zip中的所有文件
	#		GRCh38.tx2gene: 转录本和基因的对应关系文件，第一列是转录本，第二列是其所属的基因，用于后续把转录本定量转换为基因定量

	cd ${wd}
	# 列出salmon的输出文件
	find . -name quant.sf
	# ./trt_N080611/trt_N080611.salmon.count/quant.sf
	# ./trt_N061011/trt_N061011.salmon.count/quant.sf
	# ./untrt_N61311/untrt_N61311.salmon.count/quant.sf

	# 这个压缩包下载解压到本地
	zip quant.sf.zip `find . -name quant.sf`

	# 生成一个两列文件方便R导入
	# xargs接收上一步的输出，按批次提供给下游程序作为输入
	# -i: 用{}表示传递的值
	cut -f 1 sampleFile | xargs -i echo -e "{}\t{}/{}.salmon.count/quant.sf" >salmon.output
	head salmon.output
	# Samp    Samp/Samp.salmon.count/quant.sf
	# untrt_N61311    untrt_N61311/untrt_N61311.salmon.count/quant.sf
	# untrt_N052611   untrt_N052611/untrt_N052611.salmon.count/quant.sf

	# 如果没有GTF文件，可以用其他文件，只需获取转录本和基因名字对应关系就可以
	# 如果不知道对应关系，也可以把每个转录本当做一个基因进行分析
	# Trinity拼装时会生成这个文件
	# 注意修改$14, $10为对应的信息列，
	# tx2gene为一个两列文件，第一列是转录本没名字，第二列是基因名字。
	sed 's/"/\t/g' genome/GRCh38.gtf | awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) print "TXname\tGene"; if($3=="transcript") print $14, $10}' >GRCh38.tx2gene
	head GRCh38.tx2gene
	# TXname  Gene
	# ENST00000608838 ENSG00000178591
	# ENST00000382410 ENSG00000178591
	# ENST00000382398 ENSG00000125788
	# ENST00000542572 ENSG00000125788

	#tx2gene命令解析

	head genome/GRCh38.gtf

	head -n 1 genome/GRCh38.gtf | sed 's/\t/\n/g' | sed = | sed 'N;s/\n/\t/'

	sed -n '2p' genome/GRCh38.gtf | sed 's/\t/\n/g' | sed = | sed 'N;s/\n/\t/'

	head genome/GRCh38.gtf | sed 's/"/\t/g' 

	head genome/GRCh38.gtf | sed 's/"/\t/g' | tail -n 1 | tr '\t' '\n' | sed = | sed 'N;s/\n/\t/'


至此就完成了基于Salmon的所有样本基因和转录本的定量。然后下载sampleFile、GRCh38.tx2gene、salmon.output、quant.sf.zip文件到本地进行下游分析。

	# R代码，在本地windows运行
	# library("tximport")
	# library("readr")
	# salmon_file <- read.table("salmon.output", header=T,  row.names=1, sep="\t")
	# tx2gene <- read.table("genome/GRCh38.tx2gene", header=T, row.names=NULL, sep="\t")
	# txi <- tximport(salmon_file, type = "salmon",  tx2gene = tx2gene)
	# dds <- DESeqDataSetFromTximport(txi,  sample,  ~conditions)




	# 合并转录本表达量

	#cd ${wd}
	#paste `find . -name *.salmon.gene.count.tab` | \
	#  awk 'BEGIN{OFS=FS="\t" }{line=$1; \
	#    for(i=2;i<=NF;i++) if(i%2==0) {if(FNR==1) count=$i; else count=int($i+0.5); line=line"\t"count;} print line;}' \
	#  >ehbio_trans.Count_matrix.xls
	#head ehbio_trans.Count_matrix.xls

	#--------------------------------------------------------------------------------------------------
	#--------------------------------------------------------------------------------------------------



