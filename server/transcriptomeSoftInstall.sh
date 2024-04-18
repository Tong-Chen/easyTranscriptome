[TOC]

# 转录组软件安装 {#}

下面提供了4种不同平台或不同方式的转录组环境配置方法。

1. 如果自己有服务器，推荐第二种方法；
2. 如果自己windows笔记本性能比较强或有windows系统的工作站，推荐第一种放法。
3. 第3和4种方法适合想从头安装，学习整个过程。

## 在Win10安装的Ubuntu下配置转录组环境

    # Win10 安装Ubuntu系统

    # 参考链接：https://mp.weixin.qq.com/s?__biz=MzUzMjA4Njc1MA==&mid=2247491728&idx=1&sn=5529fb60cfb7fc3150368a79448d51ab&scene=21#wechat_redirect
	#  Win11还需要安装 https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi


    # 启动进入Win10的Ubuntu系统后执行下面的命令

    mkdir -p ~/transcriptome/soft
    cd ~/transcriptome/soft

    # 把soft目录加 入环境变量中，新安装的软件都软链到soft目录

    # 永久设置环境变量，用于后续持续使用
    echo 'export PATH=~/transcriptome/soft:~/miniconda3/bin:$PATH' >>~/.bash_profile

    # 更新环境变量信息
    source ~/.bash_profile
    # 临时设置环境变量，用于此次软件安装
    export PATH=~/transcriptome/soft:~/miniconda3/bin:$PATH

    # 下载conda（下载失败可以去QQ群文件Linux文件夹下载）
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f
    
    # 如果是从QQ群下载miniconda，下载后，放到c盘根目录下
    # 运行下面一行语句
    # bash /mnt/c/Miniconda3-latest-Linux-x86_64.sh -b -f 

    # 加载环境
    ~/miniconda3/condabin/conda init
    source ~/.bashrc

    # Unpack已有镜像

    # 从QQ群下载transcriptome.env.tar.gz
    # 如果是Win10+ubuntu则 放置在 C盘根目录下

    # 新建文件夹存放transcriptome环境
    mkdir -p ~/miniconda3/envs/transcriptome

    # 解压环境
    tar -vxzf /mnt/c/transcriptome.env.tar.gz -C ~/miniconda3/envs/transcriptome

    # 激活环境
    source ~/miniconda3/envs/transcriptome/bin/activate
    conda-unpack

## 在Linux/Unix服务器进行下面的操作

    mkdir -p ~/transcriptome/soft
    cd ~/transcriptome/soft

    # 把soft目录加大环境变量中，新安装的软件都软链到soft目录

    # 永久设置环境变量，用于后续持续使用
    echo 'export PATH=~/transcriptome/soft:~/miniconda3/bin:$PATH' >>~/.bash_profile

    # 更新环境变量信息
    source ~/.bash_profile
    # 临时设置环境变量，用于此次软件安装
    export PATH=~/transcriptome/soft:~/miniconda3/bin:$PATH


    # 下载conda（下载失败可以去QQ群文件Linux文件夹下载）
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f 
    # 如果是从QQ群下载miniconda，下载后，放到家目录下
    # 然后运行下面这句话
    # bash ~/Miniconda3-latest-Linux-x86_64.sh -b -f 

    # 加载环境
    ~/miniconda3/condabin/conda init
    source ~/.bashrc
    
    # Unpack已有镜像

    # 从QQ群下载transcriptome.env.tar.gz，放置在 家目录下

    # 新建文件夹存放transcriptome环境
    mkdir -p ~/miniconda3/envs/transcriptome

    # 如果是Linux服务器则 运行下面这句
    tar -vxzf ~/transcriptome.env.tar.gz -C ~/miniconda3/envs/transcriptome

    # 激活环境
    source ~/miniconda3/envs/transcriptome/bin/activate
    conda-unpack


## 利用conda从头配置Transcriptome环境

1. 配置镜像

        conda config --add channels bioconda 
        conda config --add channels defaults
        conda config --add channels conda-forge # Highest priority
                
        # Anocanda清华镜像
        # conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2/
        # conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ 
        # conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ 
        # conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
        # 后加的通道，有最高优先级
        # conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ 
        # conda config --set show_channel_urls yes

2. 安装转录组分析所需所有软件

        #注意： conda安装的samtools不能使用samtools tview, 需要自己重新编译安装
        conda create -y -n transcriptome python=3.6 r=4.1.0
        conda activate transcriptome
    
        # source ~/miniconda3/bin/activate transcriptome
        # conda安装时注意顺序
        conda install -y samtools multiqc fastqc star 
        conda install -y stringtie trimmomatic
        conda install -y rmats
        # 注意h5py的版本
        # conda install -y h5py=2.10.0 rnasamba
        conda install -y rseqc
        conda install -y salmon
        conda install -y bedtools htseq
        conda install -y rmats2sashimiplot
        conda install -y gffread
        conda install -y ucsc-bedsort ucltrsc-bedgraphtobigwig ucsc-gtftogenepred ucsc-genepredtobed

## 安装其它软件和脚本 (每种方法都必须要安装的)

    # sra-tools 

    # 如果操作系统是centos，则运行下面这句话
    cd ~/transcriptome/soft
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
    mkdir -p sratoolkit
    tar xvzf sratoolkit.current-centos_linux64.tar.gz -c sratoolkit
    ln -s `pwd`/sratoolkit/sra-tools*/bin/* ~/transcriptome/soft

    
    # 如果操作系统是ubuntu，则运行下面这句话
    cd ~/transcriptome/soft
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    mkdir -p sratoolkit
    tar xvzf sratoolkit.current-ubuntu64.tar.gz -c sratoolkit
    ln -s `pwd`/sratoolkit/sra-tools*/bin/* ~/transcriptome/soft
    

    # 下载rnasamba的模型
    # cd ~/transcriptome/soft
    # wget -c https://raw.githubusercontent.com/apcamargo/RNAsamba/master/data/full_length_weights.hdf5
    # wget -c https://raw.githubusercontent.com/apcamargo/RNAsamba/master/data/partial_length_weights.hdf5


## 安装我们自己写的脚本 和修改后的 CPC2 (每种方法都必须要安装的)

    wget https://github.com/Tong-Chen/easyTranscriptome/archive/refs/tags/v0.02.zip
    unzip v0.02.zip
    cp -rf easyTranscriptome*/* ~/transcriptome/soft/
    cd ~/transcriptome/soft/CPC2/libs/libsvm/libsvm-3.18
    make clean && make
    ln -s ~/transcriptome/soft/CPC2/bin/* ~/transcriptome/soft/
    chmod 755 ~/transcriptome/soft/*
    export PATH=$PATH:~/transcriptome/soft/

# 备选安装方案 (前面的安装都包含了所有内容，如果有失败或其他原因才使用后面的安装)

## Trinity如果要安装最新版需要自己编译

    # conda install trinity安装的是trinity-2.8.5。已经包含在了上面的环境中

    # 最新版手动编译安装
    cd ~/transcriptome/soft/
    wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.11.0/trinityrnaseq-v2.11.0.FULL.tar.gz
    tar xvzf trinityrnaseq-v2.11.0.FULL.tar.gz
    cd trinityrnaseq-v2.11.0/
    make 
    make plugins
    export TRINITY_HOME=`pwd`/trinityrnaseq-v2.11.0

## 安装gffcompare (可选，前面conda环境已经包含了)

    mkdir -p ~/transcriptome/soft/gffcom
    cd ~/transcriptome/soft/gffcom
    git clone https://github.com/gpertea/gclib
    git clone https://github.com/gpertea/gffcompare
    cd gffcompare
    make release
    ln -sf `pwd`/gffcompare ~/transcriptome/soft

    cd ~/transcriptome/soft/gffcom
    git clone https://github.com/gpertea/gffread
    cd gffread
    make
    ln -sf `pwd`/gffread ~/transcriptome/soft


## 手动安装rmats2sashimiplot (可选，前面conda环境已经包含了)

    cd ~/transcriptome/soft
    # git clone git@github.com:Xinglab/rmats2sashimiplot.git rmats2sashimiplot
    git clone https://github.com/Xinglab/rmats2sashimiplot.git rmats2sashimiplot
    cd rmats2sashimiplot
    ./2to3.sh
    python setup.py install

## 后续仅供参考（是比较老的安装方式）

    #########逐个安装##################################################################
    #
    ## SRA toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>, 根据服务器操作系统类型下载对应的二进制编码包，下载解压放到环境变量即可使用。
    ## CentOS下地址：<https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.0/sratoolkit.2.9.0-centos_linux64.tar.gz>。
    #
    #wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.0/sratoolkit.2.9.0-centos_linux64.tar.gz
    #tar xvzf sratoolkit.2.9.0-centos_linux64.tar.gz
    #ln -s `pwd`/sratoolkit.2.9.0-centos_linux64/bin/fastq-dump ~/transcriptome/soft
    ## 若运行成功，则输出为 ~/transcriptome/soft/fastq-dump
    #which fastq-dump
    #
    ## Fastqc软件安装
    ## yum install java-1.8.0-openjdk.x86_64 # 如果需要时安装java
    #cd ~/transcriptome/soft
    #wget -c http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
    #unzip fastqc_v0.11.5.zip
    #chmod 755 FastQC/fastqc
    ## 假设~/transcriptome/soft目录在环境变量中
    #ln -s `pwd`/FastQC/fastqc ~/transcriptome/soft/
    #
    #
    ## Trimmomatic安装
    #
    #cd ~/transcriptome/soft
    #wget -c http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
    #unzip Trimmomatic-0.36.zip
    #cd Trimmomatic-0.36
    #echo "java -jar \$(dirname \$(readlink -f \$0))/trimmomatic-0.36.jar \$*" >trimmomatic.sh
    #chmod 755 trimmomatic.sh
    ## 假设~/transcriptome/soft目录在环境变量中
    #ln -s `pwd`/trimmomatic.sh ~/transcriptome/soft/
    #
    ## Salmon安装
    #
    #wget -c https://github.com/COMBINE-lab/salmon/releases/download/v0.12.0-alpha/salmon-latest_linux_x86_64.tar.gz
    #tar xvzf salmon-latest_linux_x86_64.tar.gz
    #ln -s `pwd`/salmon-latest-linux_x86_64/bin/salmon ~/transcriptome/soft
    #
    ### 安装conda
    #
    #cd ~/transcriptome/soft
    #wget -c -c https://repo.anaconda.com/archive/Anaconda2-5.1.0-Linux-x86_64.sh
    #bash Anaconda2-5.1.0-Linux-x86_64.sh -b -f 
    ## 默认安装在~/anaconda2中
    #
    #
    ## cd ~/transcriptome/soft
    ## tar xvf samtools-1.6.tar
    ## cd samtools-1.6
    ## make
    ## ln -sf `pwd`/samtools ~/transcriptome/soft/
    #
    #
    #
    ### 采用STAR/HISAT2构建索引 
    #
    ##[HISAT2](http://ccb.jhu.edu/transcriptome/software/hisat2/index.shtml)安装
    #
    #cd ~/transcriptome/soft
    #wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
    #unzip hisat2-2.1.0-Linux_x86_64.zip
    #ln -s `pwd`/hisat2-2.1.0/* ~/transcriptome/soft
    #
    ## STAR
    #cd ~/transcriptome/soft
    #wget https://github.com/alexdobin/STAR/archive/2.6.0c.zip -O STAR.2.6.0c.zip
    #unzip STAR.2.6.0c.zip
    #ln -s `pwd`/STAR-2.6.0c/bin/Linux_x86_64_static/* ~/transcriptome/soft
    #
    ## RSeQC是一个python包，pip安装就可以
    #cd ~/transcriptome/soft
    #ln -s ~/anaconda2/include/lzo ~/anaconda2/include/python2.7/
    #ln -s ~/anaconda2/include/lzo/* ~/anaconda2/include/python2.7/
    #pip install -i https://pypi.tuna.tsinghua.edu.cn/simple RSeQC
    ## 下面三个脚本 geneBody_coverage2.py, read_distribution.py, RPKM_saturation.py 
    ## 都是RSeQC的一个子工具
    #
    ## HT-seq安装
    #pip install -i https://pypi.tuna.tsinghua.edu.cn/simple htseq
    #
    #
    #
    ## R包安装
    #
    ##site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
    ##options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
    #
    ### 转录本拼装
    #
    ##[StringTie](https://ccb.jhu.edu/transcriptome/software/stringtie/index.shtml?t=manual)转录本拼装
    #
    ##软件安装
    #
    ## stringtie
    #wget -c http://ccb.jhu.edu/transcriptome/software/stringtie/dl/stringtie-1.3.4d.Linux_x86_64.tar.gz
    #tar xvzf stringtie-1.3.4d.Linux_x86_64.tar.gz
    #chmod 755 stringtie-1.3.4d.Linux_x86_64/stringtie
    ## 假设~/transcriptome/soft在环境变量中
    #ln -sf `pwd`/stringtie-1.3.4d.Linux_x86_64/stringtie ~/transcriptome/soft
    #
    ## 安装prepDE.py
    ## cd ~/transcriptome/soft
    ## wget -c https://ccb.jhu.edu/transcriptome/software/stringtie/dl/prepDE.py
    ## chmod 755 prepDE.py
    ##ln -s `pwd`/prepDE.py ~/transcriptome/soft
    #
    #
    #
    ##TPM生成 RSEM安装
    ##cd ~/transcriptome/soft
    ##wget -c https://github.com/deweylab/RSEM/archive/v1.3.1.tar.gz -O RSEM.v1.3.1.zip
    ##tar xvzf RSEM.v1.3.1.zip
    ##cd RSEM-1.3.1
    ##make && make install
    ##wget -c https://raw.githubusercontent.com/miyagawa/cpanminus/master/cpanm -O ~/transcriptome/soft/cpanm
    ##cpanm Env
    #
    ## 差异剪接分析 {#alternative_splicing}
    #
    #
    ###软件安装
    #
    ## rMATS
    #pip install -i https://pypi.tuna.tsinghua.edu.cn/simple pysam
    ##yum install gsl-devel.x86_64 lapack-devel blas-devel
    #cd ~/transcriptome/soft
    #wget -c https://sourceforge.net/projects/rnaseq-mats/files/MATS/rMATS.4.0.2.tgz
    #tar xvzf rMATS.4.0.2.tgz
    #python -c "import sys; print sys.maxunicode"
    ## 输出为1114111，使用rMATS-turbo-Linux-UCS4
    #cd rMATS.4.0.2/rMATS-turbo-Linux-UCS4
    #chmod 755 rmats.py
    ## 默认使用全路径调用，如果想简化些，需要修改rmats.py源文件
    ## 第246行修改如下：
    ## root_dir = os.path.dirname(os.path.realpath(__file__))
    #ln -s `pwd`/rmats.py ~/transcriptome/soft
    #
    #
    #
    ## sashimiplot绘制
    #pip install -i https://pypi.tuna.tsinghua.edu.cn/simple rmats2sashimiplot 
    # UCSC Toolkit 下载之后，放入环境变量即可使用

    cd ~/transcriptome/soft
    wget -c ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedSort
    wget -c ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
    wget -c ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
    wget -c ftp://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed


