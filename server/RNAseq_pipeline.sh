[TOC]

    # 易转录组EasyTranscriptome
    
    # 作者 Authors: 陈同 (Chen Tong)等
    # 版本 Version: v2.0
    # 更新 Update: 2024-4-18
    # 系统要求 System requirement: Windows 10 / Mac OS 10.12+ / Ubuntu 16.04+ / CentOS
    
    # 1. 这一部分是在远程服务器的Rstudio界面利用远程服务器的计算资源进行的操作。
    # 2. 在服务器运行时，不能跳着运行，有时前一步的输出结果是后一步的输入结果。
    # 3. 文件的上传和下载演示。
    # 4. 此服务器只供练习使用，有效期一个月。
    # 5. 所有#号开头的行是注释。
    # 6. 此文为Markdown格式，放入Markdown阅读器可查看其层级结构。
    # 7. 如果代码被分割为多行（如下，行末有 \ ），建议全选运行。
    
    # 设置工作(work directory, wd)和软件/数据库(database, db)目录
    # 添加环境变量，并进入工作目录 Add environmental variables and enter work directory
    # **每次打开Rstudio必须运行下面5行 Run it**
    
    db=~/transcriptome/genome
    wd=~/transcriptome/project
    export PATH=$PATH:~/transcriptome/soft
    chmod 755 ~/transcriptome/soft/*
    cd ${wd}
    conda activate transcriptome

# 1 输入文件准备

## 1.1 测序原始数据下载 （跳过）

    # 使用NCBI提供的SRA-toolkit中的工具`fastq-dump`直接下载SRR文件，
    # 并转换为`FASTQ`格式，`--split-3`参数表示如果是双端测序就自动拆分，如果是单端不受影响。
    
    # SRA toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>, 
    # 根据服务器操作系统类型下载对应的二进制编码包，下载解压放到环境变量即可使用。
    
    # CentOS下地址：<https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.0/sratoolkit.2.9.0-centos_linux64.tar.gz>。
    
    # 如果需要下载测序数据，则取消下面语句前的#
    
    #cd ~/transcriptome/data
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
    
## 1.2 公司返回的测序结果 

    # 公司返回的测序结果，通常为一个样品对应一对fq/fastq.gz格式压缩文件
    # 文件名与metadata 中的样品名务必对应：不一致时需要手工修改，
    # 批量改名见"https://mp.weixin.qq.com/s/41mA-6q8j7EALiLDkkcyiQ"
    # -sh 表示获得每个文件的大小 (s: size)，以人可读的形式显示 (h)即加上单位如M，K，G等。
    ls -sh seq/*.fq.gz
    
    # 13M seq/trt_N052611_1.fq.gz  14M seq/trt_N61311_1.fq.gz     19M seq/untrt_N080611_1.fq.gz
    # 13M seq/trt_N052611_2.fq.gz  14M seq/trt_N61311_2.fq.gz     19M seq/untrt_N080611_2.fq.gz
    # 18M seq/trt_N061011_1.fq.gz  19M seq/untrt_N052611_1.fq.gz  15M seq/untrt_N61311_1.fq.gz
    # 18M seq/trt_N061011_2.fq.gz  19M seq/untrt_N052611_2.fq.gz  15M seq/untrt_N61311_2.fq.gz
    # 25M seq/trt_N080611_1.fq.gz  16M seq/untrt_N061011_1.fq.gz
    # 25M seq/trt_N080611_2.fq.gz  16M seq/untrt_N061011_2.fq.gz

## 1.3 测序原始数据查看和测序量评估

    # 输入文件：测序公司返回的FASTQ文件 *.fq.gz 如trt_N061011_1.fq.gz
    # 输出文件：html文件 如 trt_N061011_1_fastqc.html
    
    ## 测序质量评估 (NGS基础 - FASTQ格式解释和质量评估: https://mp.weixin.qq.com/s/tDMih7ISLJcL4F4sWBq3Vw)
    
    # 查看文件，压缩或未压缩的都可以，尤其适合查看极大文件
    # 在Rstudio界面不可用, 可以查看，但不能退出
    #less trt_N061011_1.fq.gz
    
    # zcat查看gzip压缩的文件
    # head -n 8 显示前8行文件内容
    # 前8行是代表几条序列？
    # gzip: stdout: Broken pipe: 不是错误，忽略即可 *******
    zcat seq/trt_N061011_1.fq.gz | head -n 8
    
    ### @SRR1039521.13952745/1
    ### TTCCTTCCTCCTCTCCCTCCCTCCCTCCTTTCTTTCTTCCTGTGGTTTTTTCCTCTCTTCTTC
    ### +
    ### HIJIIJHGHHIJIIIJJJJJJJJJJJJJJJJJJJJJIIJJFIDHIBGHJIHHHHHHFFFFFFE
    ### @SRR1039521.7213571/1
    ### GAGGAAGGGCAGAGGGAGCAGGGAGACTGTAGATCAGGGGCTGAATGGAGATCCGGTCCTG
    ### +
    ### <1E:E@3B:8?E?B;6FFFI@@6/0@-B<FFD4)=C=D>'-54@B@>>C>AC;>=38=6A@
    
    zcat seq/trt_N061011_2.fq.gz | head -n 8
    
    # 测序reads数计算
    
    # wc -l: 计算行数
    # bc -l: 计算器 (-l：浮点运算)
    # 为什么除以4，又除以1000000？
    echo "`zcat seq/trt_N061011_1.fq.gz | wc -l` / (4*1000000)" | bc -l
    # .464329 million
    
    # 测序碱基数计算 (awk的介绍见[常用和不太常用的awk命令](http://mp.weixin.qq.com/s/8wD14FXt7fLDo1BjJyT0ew))
    
    # awk运算
    # %取余数
    zcat seq/trt_N061011_1.fq.gz | awk '{if(FNR%4==0) base+=length}END{print base/10^9,"G";}'
    # 0.028816 G
    
    # 批量统计测序数据并汇总表
    seqkit stat seq/*.fq.gz > seq/seqkit.txt
    head seq/seqkit.txt
    
    # 批量评估，结果输出在 seq 目录下
    fastqc seq/*.fq.gz
    
    # seq/trt_N061011_1_fastqc.html 在Rstudio中打开 (View in Web browser)
    
    # multiqc 整理评估结果
    # -d .: 表示分析当前目录所有文件和文件夹，multiqc会遍历读取，分析每一个文件是不是
    #       常用工具的输出结果文件，再进行读取
    multiqc -d seq/ -o multiqc
    
    # 结果存储在multiqc/multiqc_report.html, 双击打开

## 1.4 实验设计文件生成 (metadata.txt)

    # 实验设计文件是用于将样本名和分组信息对应的文件，用途有 2：
    #   1. 利用实验设计文件中的样本名循环获得所有样品的测序reads
    #   2. 利用实验设计文件中的样本分租信息进行后续差异分析
    # 实验设计文件是TAB键分割的样品分组信息，可以用记事本、Notepad 或 Excel 等手写。
    # 如果是 Excel 写的，Excel 写完后需另存为单个TAB键（制表符）分割的文本文件，
    # 第一列是样本名字
    # 后面的列是样本分组或其它属性信息，如实验批次等，记录越全越好
    head result/metadata.txt
    # csvtk统计表行(样本数，不含表头)列数，-t设置列分隔为制表符
    csvtk -t stat result/metadata.txt
    # windows用户如果结尾有^M，运行sed命令去除，再用cat -A检查结果
    sed -i 's/\r//' result/metadata.txt
    cat -A result/metadata.txt | head -n3


# 2 基于Salmon的定量流程 {#salmon}

## 2.1 构建基因组/转录本索引 (每个物种通常只需构建一次)

### 2.1.1 基因注释文件准备

    # 工作目录：~/transcriptome/genome
    # 输入文件：
    #   Genome.gtf 人基因注释序列，从Ensembl下载
    #   Genome.fa 人基因组序列，从Ensembl下载
    # 输出文件：
    #   Genome.transcript.fa 所有转录本序列文件
    ## 获取转录本序列
    
    # 如果物种有cDNA序列提供，直接下载；但Biomart直接下载的cDNA只包含蛋白编码基因，
    # 建议自己提取
    
    # 如果是三代测序，可直接用公司提供的序列作为cDNA
    
    # 如果是无参转录组，直接用Trinity拼装的转录本作为cDNA
    
    # 如果没有则用下面的程序，通过gtf和genome文件提取
    
    # Genome.fa 人基因组序列，从Ensembl下载
    # Genome.gtf 人基因注释序列，从Ensembl下载
    # 
    # wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz -O genome.fa.gz
    # wget http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O genome.fa.gz
    # gunzip -c genome.fa.gz >Genome.fa
    # wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O gene_annotation.gtf.gz
    # gunzip -c gene_annotation.gtf.gz >Genome.gtf
    
    gunzip ${db}/Genome.fa.gz
    gunzip ${db}/Genome.gtf.gz
    gffread ${db}/Genome.gtf -g ${db}/Genome.fa -w ${db}/Genome.transcript.fa.tmp
    
    # gffread生成的fasta文件同时包含基因名字和转录本名字
    head -n 80 ${db}/Genome.transcript.fa.tmp | grep '>'
    
    # 去掉空格后面的字符串，保证cDNA文件中fasta序列的名字简洁，不然后续会出错
    cut -f 1 -d ' ' ${db}/Genome.transcript.fa.tmp >${db}/Genome.transcript.fa
    head ${db}/Genome.transcript.fa
    
    head -n 80 ${db}/Genome.transcript.fa | grep '>'

### 2.1.2 无基因组的情况下构建Salmon索引 (cDNA) 

    # 工作目录：~/transcriptome/genome
    # 输入文件：
    #   Genome.transcript.fa 所有转录本序列文件
    # 输出文件：
    #   Genome.salmon：一个文件夹，构建的salmon索引
            
    salmon index -t ${db}/Genome.transcript.fa -i ${db}/Genome_index.salmon
    
    # salmon索引目录下的内容
    # -s: 显示文件大小
    # -h: human-readable 人类可读的方式显示文件大小
    ls -sh ${db}/Genome_index.salmon
    
    # 总用量 14M
    #  20K complete_ref_lens.bin   1.9M mphf.bin              12K ref_indexing.log
    # 588K ctable.bin              7.8M pos.bin               20K reflengths.bin
    #  36K ctg_offsets.bin         4.0K pre_indexing.log     1.7M refseq.bin
    # 4.0K duplicate_clusters.tsv  424K rank.bin             844K seq.bin
    # 4.0K info.json                36K refAccumLengths.bin  4.0K versionInfo.json
    
    # 有些转录本名字不同，但序列完全一样，被鉴定出来只保留一个。
    head ${db}/Genome_index.salmon/duplicate_clusters.tsv 
    
    # RetainedTxp    DuplicateTxp
    # ENST00000486475    ENST00000629881
    # ENST00000401322    ENST00000637366
    # ENST00000401322    ENST00000401376
    # ENST00000401322    ENST00000621667

### 2.1.3 有基因组的情况下构建Salmon索引 (cDNA+genome)

    # 区分序列是来源于转录本还是其它基因组区域
    
    # 参考
    # https://github.com/COMBINE-lab/salmon
    # https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    
    # 获取所有基因组序列的名字存储于decoy中
    grep '^>' ${db}/Genome.fa | cut -d ' ' -f 1 | \
      sed 's/^>//g' >${db}/Genome.decoys.txt
    
    # 合并cDNA和基因组序列一起
    # 注意：cDNA在前，基因组在后 *********
    
    cat ${db}/Genome.transcript.fa ${db}/Genome.fa >${db}/Genome_trans_dna.fa
    ls -sh ${db}/Genome.transcript.fa ${db}/Genome.fa ${db}/Genome_trans_dna.fa
    
    # 构建索引 （更慢，结果会更准） # --no-version-check
    salmon index -t ${db}/Genome_trans_dna.fa -d ${db}/Genome.decoys.txt -i ${db}/Genome_index.salmon

    # salmon index --help

## 2.2 单个样品Salmon定量获得基因表达表 (练习用)

    # 工作目录：~/transcriptome/project/
    # 输入文件：
    #   seq/trt_N061011_1.fq.gz, seq/trt_N061011_2.fq.gz: 双端测序数据
    #   ${db}/Genome_index.salmon：基因组索引 (上一步构建获得的)
    
    # 输出文件：
    #   temp/trt_N061011/trt_N061011.salmon.count: salmon定量结果输出目录
    #   temp/trt_N061011/trt_N061011.salmon.count/quant.sf: 样品转录组定量结果，差异分析时会用到
    
    # 注意切换目录
    cd ${wd}
    
    # -p: 表示若待创建的文件夹已存在则跳过；若不存在，则创建；也可用于创建多层文件夹
    # man mkdir 可查看详细帮助
    mkdir -p result/trt_N061011
    
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
    salmon quant --gcBias -l A -1 seq/trt_N061011_1.fq.gz -2 seq/trt_N061011_2.fq.gz \
      -i ${db}/Genome_index.salmon -o result/trt_N061011/trt_N061011.salmon.count -p 1
    
    # 如果是单端数据，只需提供 -r
    # salmon quant --gcBias -l A -r seq/trt_N061011.fq.gz -i ${db}/Genome_index.salmon \
    #  -o result/trt_N061011/trt_N061011.salmon.count -p 4
    
    # 输出结果存储在 temp/trt_N061011/trt_N061011.salmon.count目录中
    # quant.sf 为转录本表达定量结果，第4列为TPM结果，第5列为reads count
    # quant.genes.sf 为基因表达定量结果
    head -n 30 result/trt_N061011/trt_N061011.salmon.count/quant.sf | tail

## 2.3 所有样本批量定量

    # 工作目录：~/transcriptome/project/
    # 输入文件：
    #   metadata.txt: TAB键分割的样品分组信息，至少2列，样本名字和样本所属分组信息
    #                 注意：样本名字与测序数据文件名前缀一致 **
    #   ${db}/Genome_index.salmon：基因组索引 (上一步构建获得的)
    # 隐式输入文件：
    #   每个样本质控后的测序数据 (对Salmon来讲，提供raw data和clean data差别不大，
    #   因为是基于K-mer的计算)
    #   metadata.txt中样本名是 trt_N061011，加上后缀_1.fq.gz可获得左端测序数据名字，
    #   加上后缀_2.fq.gz可获得右端测序数据名字。
    #   trt_N061011_1.fq.gz, trt_N061011_2.fq.gz: 双端测序数据
    # 输出文件：
    #   每个样本一个文件夹
    #     样本名/样本名.salmon.count: salmon定量结果输出目录
    #     样本名/样本名.salmon.count/quant.sf: 样品转录组定量结果，差异分析时会用到
    #       trt_N061011/trt_N061011.salmon.count: salmon定量结果输出目录
    #       trt_N061011/trt_N061011.salmon.count/quant.sf: 样品转录组定量结果，
    #                             差异分析时会用到


    # 批量定量
    # -g ${db}/Genome.gtf: 同时输出基因的定量结果，默认只有转录本定量结果
    # 注意末尾的 &，把程序放入后台，避免被取词工具的复制键 (ctrl+c)干扰。
    # ctrl+c在命令行下是打断程序运行
    # 如果代码被分割为多行（如下，行末有 \ ），建议全选运行
    
    for samp in `tail -n +2 result/metadata.txt | cut -f 1`; do \
      mkdir -p result/${samp}; \
      salmon quant --gcBias -l A -1 seq/${samp}_1.fq.gz \
      -2 seq/${samp}_2.fq.gz -i ${db}/Genome_index.salmon \
      -o result/${samp}/${samp}.salmon.count -p 2 \
      >result/${samp}/${samp}.salmon.log 2>&1; done &
      
    # 如果是单端，用下面的代码
    # for samp in `tail -n +2 result/metadata.txt | cut -f 1`; do \
    #   mkdir -p result/${samp}; \
    #   salmon quant --gcBias -l A -r seq/${samp}.fq.gz \
    #   -i ${db}/Genome_index.salmon \
    #   -o result/${samp}/${samp}.salmon.count -p 2 \
    #   >result/${samp}/${samp}.salmon.log 2>&1; done &
    
    # 命令拆解，便于理解
    # tail -n +2 result/metadata.txt
    # tail -n +2 result/metadata.txt | cut -f 1
    # for samp in `tail -n +2 result/metadata.txt | cut -f 1`; do echo ${samp}_1.fq.gz; done


## 2.4 利用multiqc整理salmon的输出日志

    # 获得reads比对率信息
    # -f覆盖之前的日志
    multiqc -f -d . -o multiqc/

## 2.5 整理Salmon输出用于后续差异基因分析

    # 工作目录：~/transcriptome/project/
    # 用于差异基因分析的3个文件：
    #    quant.sf.zip: 所有样本的quant.sf文件的集合
    #    salmon.output: 一个文本文件，用于记录每个样本和对应的quant.sf所在的路径，方便做差异基因分析时同时读取quant.sf.zip中的所有文件
    #    GRCh38.tx2gene: 转录本和基因的对应关系文件，第一列是转录本，第二列是其所属的基因，用于后续把转录本定量转换为基因定量
    
    # 列出salmon的输出文件
    find . -name quant.sf
    # ./trt_N080611/trt_N080611.salmon.count/quant.sf
    # ./trt_N061011/trt_N061011.salmon.count/quant.sf
    # ./untrt_N61311/untrt_N61311.salmon.count/quant.sf
    
    # 这个压缩包下载解压到本地
    (cd result; zip quant.sf.zip `find . -name quant.sf`)
    
    # 生成一个两列文件方便R导入
    # xargs接收上一步的输出，按批次提供给下游程序作为输入
    # -i: 用{}表示传递的值
    # cut -f 1 result/metadata.txt | xargs -i echo -e "{}\t{}/{}.salmon.count/quant.sf" >salmon.output
    awk 'BEGIN{OFS=FS="\t"}{print $1"\t"$1"/"$1".salmon.count/quant.sf"}' result/metadata.txt >result/salmon.output
    head result/salmon.output
    # Samp    Samp/Samp.salmon.count/quant.sf
    # untrt_N61311    untrt_N61311/untrt_N61311.salmon.count/quant.sf
    # untrt_N052611   untrt_N052611/untrt_N052611.salmon.count/quant.sf
    
    # 如果没有GTF文件，可以用其他文件，只需获取转录本和基因名字对应关系就可以
    # 如果不知道对应关系，也可以把每个转录本当做一个基因进行分析
    # Trinity拼装时会生成这个文件
    # 注意修改$14, $10为对应的信息列，
    # tx2gene为一个两列文件，第一列是转录本没名字，第二列是基因名字。
    sed 's/"/\t/g' ${db}/Genome.gtf | awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) print "TXname\tGene"; if($3=="transcript") print $14, $10}' >result/GRCh38.tx2gene
    head result/GRCh38.tx2gene
    # TXname  Gene
    # ENST00000608838 ENSG00000178591
    # ENST00000382410 ENSG00000178591
    # ENST00000382398 ENSG00000125788
    # ENST00000542572 ENSG00000125788
    
    # For mouse GTF
    # sed 's/"/\t/g' ${db}/Genome.gtf | awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) print "TXname\tGene"; if($3=="transcript") print $14, $10"_"$18}' | head


    #tx2gene命令解析
    
    head ${db}/Genome.gtf
    
    head -n 1 ${db}/Genome.gtf | sed 's/\t/\n/g' | sed = | sed 'N;s/\n/\t/'
    
    sed -n '2p' ${db}/Genome.gtf | sed 's/\t/\n/g' | sed = | sed 'N;s/\n/\t/'
    
    head ${db}/Genome.gtf | sed 's/"/\t/g' 
    
    head ${db}/Genome.gtf | sed 's/"/\t/g' | tail -n 1 | tr '\t' '\n' | sed = | sed 'N;s/\n/\t/'
    
    sed -n '2p' ${db}/Genome.gtf | sed 's/"/\t/g' | tail -n 1 | tr '\t' '\n' | awk 'BEGIN{OFS=FS="\t"}{print FNR, $0}'


至此就完成了基于Salmon的所有样本基因和转录本的定量。然后下载metadata.txt、GRCh38.tx2gene、salmon.output、quant.sf.zip文件到本地进行下游分析。


# 3  基于STAR进行转录组比对和定量的流程 {#STAR}

## 3.1 STAR 构建基因组索引 （每个物种通常只需构建 1 次索引）

### 3.1.1  基因注释文件准备和格式转换

    ## 基因注释文件准备 gtf转bed12
    
    # gtf2bed12.sh所做的操作就是下面两句话
    # wget -c ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedSort
    # wget -c ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
    # wget -c ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
    # wget -c ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
    # gtfToGenePred -ignoreGroupsWithoutExons GRCh38.gtf GRCh38.gtf.50505050.pred
    # genePredToBed GRCh38.gtf.50505050.pred GRCh38.gtf.bed12
    # 删除临时文件
    # /bin/rm -f GRCh38.gtf.50505050.pred
    
    gtf2bed12.sh -f ${db}/Genome.gtf
    # 选择长度适中的转录本用于后续评估
    awk '$3-$2>1000 && $3-$2<2000' ${db}/Genome.gtf.bed12 >${db}/Genome.model.gtf.bed12
    head ${db}/Genome.model.gtf.bed12
    
    # 生成igv截图脚本
    cut -f 1-6 ${db}/Genome.gtf.bed12 >result/igv.capture.bed
    bedtools igv -i result/igv.capture.bed -path "D:" >result/igv.batch.script
    
    # # cp ../bak/${db}/Genome.chromsize .
    # 下载chrome size文件，第一列为染色体名字，第二列为染色体大小
    # 也可以自己写程序计算
    # 适用于任何基因组序列，先判断是不是读到的第一条染色体，如果不是就打印；如果是只存储
    # 记得输出最后一条染色体的长度信息
    # awk 'BEGIN{OFS="\t"}{if($0~/>/) {if(size>0) print chrname, size; size=0; chrname=$0; sub(">","", chrname);} else size+=length;}END{print chrname,size}' ${db}/Genome.fa
    # mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg38.chromInfo" | tail -n +2 >${db}/Genome.chromsize
    # head ${db}/Genome.chromsize
    
### 3.1.2 STAR构建基因组索引 (基因组索引构建的目的是为了方便快速比对)
    
    
    # 需要GTF文件中有第三列的区域类型中有exon，第9列的区域属性信息中有gene_id\transcript_id
    
    mkdir -p ${db}/Genome_index.STAR
    # --runThreadN 2: 指定使用2个线程
    # --sjdbOverhang 100: 默认
    
    # 如果代码被分割为多行（如下，行末有 \ ），建议全选运行
    STAR --runMode genomeGenerate --runThreadN 2 --genomeDir ${db}/Genome_index.STAR \
        --genomeFastaFiles ${db}/Genome.fa --sjdbGTFfile ${db}/Genome.gtf &
    
    ## 等待 Terminal 右上角的红点消失后再继续运行 ****** 
    
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
    
    ls -sh ${db}/Genome_index.STAR
    
    # 总用量 2.1G
    # 4.0K chrLength.txt      368K exonInfo.tab          1.5G SAindex
    # 4.0K chrNameLength.txt   24K geneInfo.tab          204K sjdbInfo.txt
    # 4.0K chrName.txt         64M Genome                204K sjdbList.fromGTF.out.tab
    # 4.0K chrStart.txt       4.0K genomeParameters.txt  204K sjdbList.out.tab
    # 732K exonGeTrInfo.tab   516M SA                    224K transcriptInfo.tab
    
    cp ${db}/Genome_index.STAR/chrNameLength.txt ${db}/Genome.chromsize
    
    head ${db}/Genome.chromsize
    
    # STAR解析后的基因数
    head -n 1 ${db}/Genome_index.STAR/geneInfo.tab
    
    # 原始GTF的基因数
    grep -cP '\tgene\t' ${db}/Genome.gtf
    
## 3.2  基于 STAR 的 reads 比对和定量 (单个样品比对，只用做练习)
    
### 3.2.1  基于 STAR 的 reads 比对和定量命令
    
  
    # mkdir新建目录
    # -p: 表示目录不存在则新建；存在则不做认可事情 
    mkdir -p result/trt_N061011
    # --runThreadN 4: 使用4个线程
    # --readFilesIn: 输入文件，左端和右端
    # --readFilesCommand zcat：gzip压缩，指定解压方式
    # --genomeDir：基因组索引目录的位置
    # -S: 指定输出文件
    # --genomeLoad LoadAndKeep : 共享内存
    # 其他参数自己对着star的帮助手册看
    
    # STAR比对单个样品 
    # # 动物一般写 1000000，植物一般写5000
    # 对应参数 alignIntronMax 和 alignMatesGapMax
    # max_intron_size=1000000
    STAR --runMode alignReads --runThreadN 4 \
        --readFilesIn seq/trt_N061011_1.fq.gz seq/trt_N061011_2.fq.gz \
        --readFilesCommand zcat --genomeDir ${db}/Genome_index.STAR \
        --outFileNamePrefix result/trt_N061011/trt_N061011.  \
        --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
        --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
        --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
        --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
        --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts
    
    # 如果是单端序列，用下面的命令
    #  STAR --runMode alignReads --runThreadN 4 \
    #     --readFilesIn seq/trt_N061011.fq.gz \
    #     --readFilesCommand zcat --genomeDir ${db}/Genome_index.STAR \
    #     --outFileNamePrefix result/trt_N061011/trt_N061011.  \
    #     --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    #     --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    #     --alignIntronMin 20 --alignIntronMax 1000000 \
    #     --alignMatesGapMax 1000000 \
    #     --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
    #     --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
    #     --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
    #     --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
    #     --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts
    
    # Aug 03 11:44:27 ..... started STAR run
    # Aug 03 11:44:27 ..... loading genome
    # Aug 03 11:44:30 ..... started mapping
    # Aug 03 11:44:48 ..... finished successfully
    
    
### 3.2.2  基于 STAR 的 reads 比对和 定量结果查看
    
    # 输出文件
    
    ls -sh result/trt_N061011/*
    
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
    
### 3.2.3  基于 STAR 的 reads 比对和定量结果格式转换

    #筛选reads，按坐标排序、索引BAM文件供下游使用, 也可导入IGV查看reads比对情况、是否有变异位点等。
    
    # samtools具体参数解释见 samtools -?
    # -@ 4: 4个线程
    mkdir -p tmp
    samtools sort -@ 4 -T /tmp/trt_N061011_$(whoami) \
        -o result/trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
        result/trt_N061011/trt_N061011.Aligned.out.bam
    samtools index result/trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam
    # 为什么要按坐标排序？
    # 为什么要建索引？
    # 就可以导入IGV中查看reads的比对情况了
    
    # BigWig峰图文件生成，导入IGV或UCSC genomebrowser获取表达丰度图。
    
    # Wig里面有什么？
    # 为什么要生成BigWig？
    # 是否需要标准化？
    
    STAR --runMode inputAlignmentsFromBAM \
        --inputBAMfile result/trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
        --outWigType bedGraph --outFileNamePrefix result/trt_N061011/trt_N061011. \
        --outWigNorm RPM --outWigStrand Unstranded
    bedSort result/trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bg \
        result/trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bg
    bedGraphToBigWig result/trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bg \
        ${db}/Genome_index.STAR/chrNameLength.txt \
        result/trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bw
    
    # bam2wig.py是RseQC中的软件，可以使用`which bam2wig.py`查看命令的路径
    # -s: 是chrom size文件，两列文件，前面又介绍
    # -t: 归一化因子，所有样品都归一化到统一测序深度，方便比较
    # bam2wig.py -i result/trt_N061011/trt_N061011.sortP.bam -s ${db}/Genome.chromsize \
    #     -o result/trt_N061011/trt_N061011 -t 1000000000 -q 0
    # ls -ltr result/trt_N061011/trt_N061011*
    

### 3.2.4  基于 STAR 的 reads 比对和定量结果评估

    ## MultiQC查看STAR比对结果
    
    # multiqc -f -d . -o multiqc
    
    ### 比对质量评估
    #
    #### Reads在基因上的分布评估
    #
    ## 程序运行结束后，默认生成折线图
    ## 如果样品多，也可以生成理论课件中的热图
    # 出现 WARNING: ignoring environment value of R_HOME *** 不用理会
    geneBody_coverage2.py -i \
        result/trt_N061011/trt_N061011.Signal.UniqueMultiple.str1.out.bw \
        -r ${db}/Genome.model.gtf.bed12 -o result/trt_N061011/trt_N061011.geneBody_coverage
    #
    #### Reads比对到基因组标志性区域的分布
    #
    ## 请完成堆积柱状图的绘制
    ## 或使用www.ehbio.com/ImageGP
    read_distribution.py -i result/trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
        -r ${db}/Genome.gtf.bed12 >result/trt_N061011/trt_N061011.read_distrib.xls
    #
    cat result/trt_N061011/trt_N061011.read_distrib.xls
    #
    #
    #### 测序饱和度评估
    #
    ## -s: 采样频率，0-100之间的整数，类似于步长
    ## -q: 过滤低质量比对
    RPKM_saturation.py -i \
        result/trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
        -r ${db}/Genome.gtf.bed12 -s 10 -q 0 -o result/trt_N061011/trt_N061011.RPKM_saturation
    #
    #ls -ltr trt_N061011
    #---------------------------------------------------------
    
    # multiqc -f -d . -o multiqc
    
## 3.3  基于 STAR 的 批量比对和定量

### 3.3.1  基于 STAR 的 批量比对和定量

    # 测试是否为TAB键分割
    cat -A result/metadata.txt
    # Samp^Iconditions$
    # untrt_N61311^Iuntrt$
    # untrt_N052611^Iuntrt$
    # untrt_N080611^Iuntrt$
    # untrt_N061011^Iuntrt$
    # trt_N61311^Itrt$
    # trt_N052611^Itrt$
    # trt_N080611^Itrt$
    # trt_N061011^Itrt$
    
    # 动物一般写 1000000，植物一般写5000
    # max_intron_size=1000000
    
    for i in `tail -n +2 result/metadata.txt | cut -f 1`; do 
        mkdir -p result/${i}
        STAR --runMode alignReads --runThreadN 4 \
            --readFilesIn seq/${i}_1.fq.gz seq/${i}_2.fq.gz \
            --readFilesCommand zcat --genomeDir ${db}/Genome_index.STAR \
            --outFileNamePrefix result/${i}/${i}. --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
            --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
            --alignIntronMin 20 --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
            --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
            --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
            --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
            --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts
        samtools sort -@ 10 -T ${i}.tmp \
            -o result/${i}/${i}.Aligned.sortedByCoord.out.bam \
            result/${i}/${i}.Aligned.out.bam
        samtools index result/${i}/${i}.Aligned.sortedByCoord.out.bam
        STAR --runMode inputAlignmentsFromBAM \
            --inputBAMfile result/${i}/${i}.Aligned.sortedByCoord.out.bam \
            --outWigType bedGraph --outFileNamePrefix result/${i}/${i}. \
            --outWigNorm RPM --outWigStrand Unstranded
        bedSort result/${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
            result/${i}/${i}.Signal.UniqueMultiple.str1.out.bg
        bedGraphToBigWig result/${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
            ${db}/Genome_index.STAR/chrNameLength.txt \
            result/${i}/${i}.Signal.UniqueMultiple.str1.out.bw
    done &
    
    # # 单端测序的代码
    # for i in `tail -n +2 result/metadata.txt | cut -f 1`; do 
    #     mkdir -p result/${i}
    #     STAR --runMode alignReads --runThreadN 4 \
    #         --readFilesIn seq/${i}.fq.gz \
    #         --readFilesCommand zcat --genomeDir ${db}/Genome_index.STAR \
    #         --outFileNamePrefix result/${i}/${i}. --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    #         --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    #         --alignIntronMin 20 --alignIntronMax 1000000 \
    #         --alignMatesGapMax 1000000 \
    #         --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
    #         --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
    #         --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
    #         --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
    #         --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts
    #     samtools sort -@ 10 -T ${i}.tmp \
    #         -o result/${i}/${i}.Aligned.sortedByCoord.out.bam \
    #         result/${i}/${i}.Aligned.out.bam
    #     samtools index result/${i}/${i}.Aligned.sortedByCoord.out.bam
    #     STAR --runMode inputAlignmentsFromBAM \
    #         --inputBAMfile result/${i}/${i}.Aligned.sortedByCoord.out.bam \
    #         --outWigType bedGraph --outFileNamePrefix result/${i}/${i}. \
    #         --outWigNorm RPM --outWigStrand Unstranded
    #     bedSort result/${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
    #         result/${i}/${i}.Signal.UniqueMultiple.str1.out.bg
    #     bedGraphToBigWig result/${i}/${i}.Signal.UniqueMultiple.str1.out.bg \
    #         ${db}/Genome_index.STAR/chrNameLength.txt \
    #         result/${i}/${i}.Signal.UniqueMultiple.str1.out.bw
    # done &    
    
### 3.3.2  STAR比对和定量结果的批量评估
    
    ## 批量评估
    
    # cd ~/transcriptome/data
    # for i in `tail -n +2 sampleFile | cut -f 1`; do 
    #   /anaconda2/bin/geneBody_coverage2.py -i \
    #     ${i}/${i}.Signal.UniqueMultiple.str1.out.bw \
    #     -r genome/GRCh38.model.gtf.bed12 -o ${i}/${i}.geneBody_coverage
    #   /anaconda2/bin/read_distribution.py -i ${i}/${i}.Aligned.sortedByCoord.out.bam \
    #     -r ${db}/Genome.gtf.bed12 >${i}/${i}.read_distrib.xls
    #   /anaconda2/bin/RPKM_saturation.py -i \
    #     ${i}/${i}.Aligned.sortedByCoord.out.bam \
    #     -r ${db}/Genome.gtf.bed12 -s 10 -q 0 -o ${i}/${i}.RPKM_saturation
    # done
    
    ## multiqc整理软件运行结果
    
    multiqc -f -d . -o multiqc
    
## 3.4  合并STAR 比对结果获得 Count 矩阵
    
    # 基因reads count增加样品信息，方便后续合并        
    # sed '5 i\Gene\ttrt_N061011\ttrt_N061011\ttrt_N061011' trt_N061011/trt_N061011.ReadsPerGene.out.tab trt_N061011/trt_N061011.ReadsPerGene.out.tab.ehbio
    
    for i in `tail -n +2 result/metadata.txt | cut -f 1`; do 
        sed "5 i\Gene\t${i}\t${i}\t${i}" result/${i}/${i}.ReadsPerGene.out.tab \
          >result/${i}/${i}.ReadsPerGene.out.tab.ehbio
    done
    
    
    # Linux 命令合并
    paste `find . -name *.ReadsPerGene.out.tab.ehbio` | tail -n +5 | \
        awk 'BEGIN{OFS=FS="\t" }{line=$1; \
        for(i=2;i<=NF;i++) if(i%2==0 && i%4!=0) line=line"\t"$i; print line;}' \
        >result/ehbio_trans.Count_matrix.xls
    head result/ehbio_trans.Count_matrix.xls
    
    
    
# 4 鉴定新基因或转录本和可变剪接 {#stringtie}

## 4.1 转录本拼装 基于 STAR 的比对结果
    
    使用STAR比对的结果拼装时，一定要加比对参数`--outSAMattrIHstart 0 --outSAMstrandField intronMotif`，不然出来的都是单外显子转录本。具体比对见上一步。
    
## 4.2 单样本拼装 (练习用)
    
    # # -G: 指定reference GTF
    # -p 1: 使用一个线程，多核处理器可调大
    # -f 0.01 : 允许的最小isoform比例，默认0.01
    
    stringtie result/trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
        -G ${db}/Genome.gtf -l trt_N061011 -o result/trt_N061011/trt_N061011.stringtie_first.gtf \
        -f 0.01 -p 2
        
    head result/trt_N061011/trt_N061011.stringtie_first.gtf
        
## 4.3  多样本批量拼装

    # 多样本循环拼装
    for i in `tail -n +2 result/metadata.txt | cut -f 1`; do 
        stringtie result/${i}/${i}.Aligned.sortedByCoord.out.bam \
                -G ${db}/Genome.gtf -l ${i} -o result/${i}/${i}.stringtie_first.gtf -f 0.01 -p 1
    done &
    
## 4.4 多样本拼装结果合并、比较，构建参考基因集
    
    # 转录本合并
    # 获取所有拼装好的gtf文件
    find . -name *.stringtie_first.gtf >result/mergeList.txt
    
    head result/mergeList.txt
    
    # -G: 指定reference GTF
    # -l: 输出结果中转录本的名字前缀
    # -o: 输出文件
    # mergeList.txt：单个样品GTF列表，每个一行
    
    stringtie --merge -G ${db}/Genome.gtf -l ehbio_trans -o result/ehbio_trans.gtf result/mergeList.txt
    
    
    #新拼装转录本与原基因组注释转录本比较 (可用来筛选新转录本)
    
    # -R: 只考虑与拼装的转录本有重叠的注释转录本
    # -r: reference gtf
    # -o: 输出前缀
    
    gffcompare -R -r ${db}/Genome.gtf -o result/assembeCompare2Ref result/ehbio_trans.gtf
    # 会输出一个assembeCompare2Ref.annotated.gtf，用于后续的定量
    
    head result/assembeCompare2Ref.annotated.gtf
    
    统计不同种类转录本的数目
    
    cut -f 3 result/assembeCompare2Ref.ehbio_trans.gtf.tmap | tail -n +2 | sort | uniq -c
    
## 4.5  新参考基因集定量 (HTSeq)

    # 转录本定量
    
    ## 基于HTSEQ
    
    for i in `tail -n +2 result/metadata.txt | cut -f 1`; do 
        htseq-count -f bam -r pos -a 10 -t exon -s no -i gene_id \
                -m union result/${i}/${i}.Aligned.sortedByCoord.out.bam \
                result/assembeCompare2Ref.annotated.gtf >result/${i}/${i}.readsCount
        grep -v '^__' result/${i}/${i}.readsCount | sed "1 iGene\t${i}"  \
                >result/${i}/${i}.readsCount2
    done &

    # 合并表达文件和差异分析同使用基因组注释文件计算差异基因一致。
    
    paste `find . -name *.readsCount2` | \
        awk 'BEGIN{OFS=FS="\t" }{line=$1; \
        for(i=2;i<=NF;i++) if(i%2==0) line=line"\t"$i; print line;}' \
        >result/ehbio_trans.Count_matrix2.xls
    head result/ehbio_trans.Count_matrix2.xls
    
    tail result/ehbio_trans.Count_matrix2.xls
    
## 4.6 新参考基因集定量 (Salmon)
    
    ## 自己完成
    
    
# 5 差异剪接分析 {#alternative_splicing}
    
## 5.1 rMATS差异剪接分析

    # -b1, -b2参数的值是一个文件，文件内是一个样品多个重复的bam文件
    # 所有bam文件写在一行，用逗号隔开
    # 下面例子中考虑到每行能展示的宽度有限，所以文件内容做了换行处理，实际是一行
    
    /bin/rm -f result/*.bam.txt
    awk 'BEGIN{OFS=FS="\t"}{if(FNR>1) a[$2]=a[$2]==""?"./result/"$1"/"$1".Aligned.sortedByCoord.out.bam":a[$2]",./result/"$1"/"$1".Aligned.sortedByCoord.out.bam"}END{for(i in a) {file_name="result/"i".bam.txt"; print a[i] >file_name;}}' result/metadata.txt
    
    head result/*.bam.txt
    
    # --gtf: 指定基因注释文件，可以是自己下载的注释，也可以是stringtie拼装完merge后的注释
    # --od: 输出目录 (output dir)
    # -t: paired-end or single-end
    # --libtype: 链特异性类型
    # -nthread: 多线程计算; --tstat: 统计模型多线程计算
    # --cstat: 差异剪接阈值，默认0.0001, 表示有0.01%的差异，取值在0-1之间
    
    mkdir -p result/rmats_trt_untrt 
    /bin/rm -rf result/rmats_trt_untrt 
    /bin/rm -rf rmats_tmp
    ## ****用于自己的数据时，result/trt.bam.txt 和 result/untrt.bam.txt需要按实际情况修改***
    rmats.py --b1 result/trt.bam.txt --b2 result/untrt.bam.txt --gtf ${db}/Genome.gtf \
          --od result/rmats_trt_untrt -t paired --libType fr-unstranded \
          --readLength 63 --nthread 2 --tstat 2 --cstat 0.0001 --tmp rmats_tmp
    # assembeCompare2Ref.annotated.gtf
    
    
    ls result/rmats_trt_untrt/
    #筛选SE差异显著的剪接位点
    
    awk 'FNR==1 || $20<0.2' result/rmats_trt_untrt/SE.MATS.JC.txt >result/rmats_trt_untrt/SE.MATS.JC.sig.txt
    
    head result/rmats_trt_untrt/SE.MATS.JC.sig.txt
    
    #筛选MXE差异显著的剪接位点
    
    awk 'FNR==1 || $20<0.2' result/rmats_trt_untrt/MXE.MATS.JC.txt >result/rmats_trt_untrt/MXE.MATS.JC.sig.txt
    
    head result/rmats_trt_untrt/MXE.MATS.JC.sig.txt

## 5.2 sashimiplot可视化差异剪接结果
    
    # --b1, --b2 同rMATS，只是直接跟文件内容而不是文件名
    # -t: 想要绘制的类型
    # -e: rMATS的输出，一般选择差异显著的绘制
    # --exon_s,  --intron_s: 缩放外显子或内含子，默认无缩放；一般用于内含子太大时，把内含子
    # 相对缩小一点，图会更好看一些
    # --l1, --l2，对应于--b1, --b2的样品的组名
    # -o 指定输出目录
    
      # group1 和 group2 根据实际修改 
      group1="trt"
      group2="untrt"

      group1_bam=$(awk -v group1=${group1} 'BEGIN{OFS=FS="\t"}{if($2==group1) a=a==""?"./result/"$1"/"$1".Aligned.sortedByCoord.out.bam":a",./result/"$1"/"$1".Aligned.sortedByCoord.out.bam"}END{print a;}' result/metadata.txt)
      group2_bam=$(awk -v group2=${group2} 'BEGIN{OFS=FS="\t"}{if($2==group2) a=a==""?"./result/"$1"/"$1".Aligned.sortedByCoord.out.bam":a",./result/"$1"/"$1".Aligned.sortedByCoord.out.bam"}END{print a;}' result/metadata.txt)

    rmats2sashimiplot --b1 ${group1_bam} --b2 ${group2_bam} \
        -t SE -e result/rmats_trt_untrt/SE.MATS.JC.sig.txt --exon_s 1 --intron_s 5 \
            --l1 trt --l2 untrt -o result/rmats_trt_untrt/SE.MATS.JS.sig.sashimiplot

    
# 6  无参转录组分析 (Trinity) 
    
    # 列出所有fastq文件
    ls seq/*.fq.gz
    
    # 合并左端序列和右端序列
    
    # cat *_1.fq.gz >Trinity_source_1.fq.gz
    # cat *_2.fq.gz >Trinity_source_2.fq.gz
    
    cat `cut -f 1 result/metadata.txt | tail -n +2 | sed -e 's#^#seq/#' -e 's/$/_1.fq.gz/'` >seq/Trinity_source_1.fq.gz
    cat `cut -f 1 result/metadata.txt | tail -n +2 | sed -e 's#^#seq/#' -e 's/$/_2.fq.gz/'` >seq/Trinity_source_2.fq.gz
    
    ls -sh seq/Trinity_source_1.fq.gz  seq/Trinity_source_2.fq.gz
    
    # Trinity拼装
    
    Trinity --full_cleanup --seqType fq --max_memory 10G --CPU 2 --output result/trinity.tmp \
        --left seq/Trinity_source_1.fq.gz --right seq/Trinity_source_2.fq.gz
    
# 7  新非编码 RNA 鉴定和 ceRNA 分析
    
## 7.1 提取新旧转录本及其序列

    # 前面获取了拼装完成的转录本
    
    head result/assembeCompare2Ref.annotated.gtf
    
    # 统计不同种类转录本的数目
    
    cut -f 3 result/assembeCompare2Ref.ehbio_trans.gtf.tmap | tail -n +2 | sort | uniq -c
    
    # 获取所有转录本的序列
    mkdir -p result2
    gffread result/assembeCompare2Ref.annotated.gtf -g ${db}/Genome.fa -w result2/Stringtie.transcript.fa.tmp
    
    # gffread生成的fasta文件同时包含基因名字和转录本名字
    head -n 300 result2/Stringtie.transcript.fa.tmp | grep '>'
    
    # 去掉空格后面的字符串，保证cDNA文件中fasta序列的名字简洁，不然后续会出错
    cut -f 1 -d ' ' result2/Stringtie.transcript.fa.tmp >result2/Stringtie.transcript.fa
    head result2/Stringtie.transcript.fa
    
    
## 7.2 定量新旧转录本
    
    salmon index -t result2/Stringtie.transcript.fa -i result2/Stringtie.transcript_index.salmon
    
    for samp in `tail -n +2 result/metadata.txt | cut -f 1`; do \
        salmon quant --gcBias -l A -1 seq/${samp}_1.fq.gz -2 seq/${samp}_2.fq.gz \
        -i result2/Stringtie.transcript_index.salmon -o result2/${samp}/${samp}.salmon.count -p 2 \
        -g result/assembeCompare2Ref.annotated.gtf >result2/${samp}.salmon.log 2>&1; done
    
    for samp in `tail -n +2 result/metadata.txt | cut -f 1`; do \
        cut -f 1,4 result2/${samp}/${samp}.salmon.count/quant.sf | \
        sed "1 s/TPM/${samp}/" >result2/${samp}/${samp}.salmon.count/quant2.sf; done
    
    paste `find . -name quant2.sf` | \
        awk 'BEGIN{OFS=FS="\t" }{line=$1; \
        for(i=2;i<=NF;i++) if(i%2==0) line=line"\t"$i; print line;}' \
        >result2/ehbio_trans.TPM.xls
    head result2/ehbio_trans.TPM.xls
    
    # 生成转录本和Gene symbol的对应关系
    
    head result/assembeCompare2Ref.annotated.gtf
    
    awk 'BEGIN{OFS=FS="\t"; read_one_line=0}{if($3=="transcript" && read_one_line==0) \
        {print $0; read_one_line=1}}' result/assembeCompare2Ref.annotated.gtf | \
        head -n 1 | sed 's/"/\t/g' | tr '\t' '\n' | sed = | sed 'N;s/\n/\t/'
        
    awk '$3=="transcript"' result/assembeCompare2Ref.annotated.gtf | \
            sed 's/"/\t/g' | awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) print "txName\tSymbol"; \
            if($14=="") $14=$10; print $10,$14}' >result2/assembeCompare2Ref.idmap.txt
    head result2/assembeCompare2Ref.idmap.txt

## 7.3 获取新转录本的序列
    
    sed 's/"/\t/g' result/assembeCompare2Ref.annotated.gtf | grep -vP 'class_code \t=\t' | \
          awk 'BEGIN{OFS=FS="\t"}{if($3=="transcript") print $10, $22}' \
          >result2/assembeCompare2Ref.annotated.new.transcript.id
    
    extractFastaByName.py -i result2/Stringtie.transcript.fa \
          -n result2/assembeCompare2Ref.annotated.new.transcript.id \
          >result2/New_transcript.fa
    
    grep '>' result2/New_transcript.fa | head
    
## 7.4 预测 lncRNA

    # 编码能力预测
    
    # rnasamba_model=~/transcriptome/soft/full_length_weights.hdf5 
    # rnasamba classify -p result2/New_transcript.translated.protein.fa \
        #  result2/New_transcript.codingstatus.txt result2/New_transcript.fa \
        #  ${rnasamba_model}
        # Change the following line in CPC2.py as
        # abspath to realpath
        # script_dir,filename = os.path.split(os.path.realpath(sys.argv[0]))
      CPC2.py -i result2/New_transcript.fa -o result2/New_transcript.codingstatus
    
    # 查看预测结果
    head result2/New_transcript.codingstatus.txt
    
    # Extract non-coding RNA
    grep 'noncoding' result2/New_transcript.codingstatus.txt >result2/New_noncoding_transcriptID.txt
    head result2/New_noncoding_transcriptID.txt
    wc -l result2/New_noncoding_transcriptID.txt
    
    extractFastaByName.py -i result2/New_transcript.fa \
          -n result2/New_noncoding_transcriptID.txt \
          >result2/New_noncoding_transcript.fa
    
    # grep: 写错误: 断开的管道 Write failed: Broken pipe
    # 不是错，忽略就好
    grep '>' result2/New_noncoding_transcript.fa | head
    
    grep -c '>' result2/New_noncoding_transcript.fa
    
    head result2/New_noncoding_transcript.fa | cut -c 1-60 
    
    # 查看翻译出的蛋白序列
    #head result2/New_transcript.translated.protein.fa
    
## 7.5 miRNA 靶点预测

    # miRNA序列下载
    
    # wget -c ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
    # gunzip -c mature.fa.gz >mature.fa
    # grep -A 1 '^>hsa-' mature.fa | sed '/--/d' | cut -f 1 -d ' ' | sed "s/hsa-//" | sed 's/U/T/g' >mirna.mature.fa
    
    head ${db}/mirna.mature.fa
    
    # miRNA-noncodingRNA靶基因预测
    
    miranda ${db}/mirna.mature.fa result2/New_noncoding_transcript.fa -sc 155 -en -20 | \
          grep -B 16 '^>>' >result2/lncRNA.miranda.predict.target
    
    grep '>>' result2/lncRNA.miranda.predict.target | cut -f 1,2 | sed 's/>>//' | \
          awk 'BEGIN{OFS=FS="\t"}{a[$1"\t"$2]+=1;}END{for(i in a) print i,a[i],"lncRNA"}' \
          >result2/lncRNA.miranda.predict.target.id
    
    head result2/lncRNA.miranda.predict.target.id
    
    # miRNA coding RNA 3UTR 靶基因预测
    
    miranda ${db}/mirna.mature.fa ${db}/UTR.fa -sc 155 -en -20 | \
          grep -B 16 '^>>' >result2/proteinCoding.miranda.predict.target
    
    grep '>>' result2/proteinCoding.miranda.predict.target | cut -f 1,2 | sed 's/>>//' | \
          awk 'BEGIN{OFS=FS="\t"}{a[$1"\t"$2]+=1;}END{for(i in a) print i,a[i],"Coding"}' >\
          result2/proteinCoding.miranda.predict.target.id
    
## 7.6 miRNA sponge 筛选
    
    awk 'BEGIN{OFS=FS="\t";}ARGIND==1{lncTar[$1]=lncTar[$1]==""?$0:lncTar[$1]"\n"$0;}\
          ARGIND==2{if(lncTar[$1]!="") {print $0; if(output[$1]=="") {print lncTar[$1]; output[$1]=1;}}}' \
          result2/lncRNA.miranda.predict.target.id result2/proteinCoding.miranda.predict.target.id \
          >result2/miRNA_sponge.list.id
    
    head result2/miRNA_sponge.list.id
    
    # Correlation 
    awk 'BGEIN{OFS=FS="\t"}ARGIND==1{a[$2]=1;}ARGIND==2{if(FNR==1 || a[$1]==1) print $0}' \
          result2/miRNA_sponge.list.id result2/ehbio_trans.TPM.xls >result2/miRNA_sponge.list.expr
    
    head result2/miRNA_sponge.list.expr
    
    pearsonCorrelation2file.py -i result2/miRNA_sponge.list.expr \
      -j result2/miRNA_sponge.list.expr \
      --positive-cor-only -o result2/miRNA_sponge.list.cor
    
    head result2/miRNA_sponge.list.cor.scipy.pearson.xls
    
    grep ENST00000246190 result2/miRNA_sponge.list.cor.scipy.pearson.xls
    
    selectSponge.py -i result2/miRNA_sponge.list.id \
          -c result2/miRNA_sponge.list.cor.scipy.pearson.xls -o result2/miRNA_sponge.pos_cor
    
    head result2/miRNA_sponge.pos_cor.table.txt
    
    head result2/miRNA_sponge.pos_cor.network.txt

