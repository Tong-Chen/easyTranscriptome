# 鉴定新基因或转录本和可变剪接 {#stringtie}

## 服务器登录

## IP: 192.168.1.107 端口 22

## 用户名: 姓名汉语拼音

## 密码: yishengxin

# 设置工作目录
wd=/mnt/d/train/serverData/data
  
# 激活环境
conda activate transcriptome

cd ${wd}



# 前面获取了拼装完成的转录本

head assembeCompare2Ref.annotated.gtf

# 统计不同种类转录本的数目

cut -f 3 assembeCompare2Ref.ehbio_trans.gtf.tmap | tail -n +2 | sort | uniq -c



# 获取所有转录本的序列

gffread assembeCompare2Ref.annotated.gtf -g genome/GRCh38.fa -w GRCh38.transcript.fa.tmp

# gffread生成的fasta文件同时包含基因名字和转录本名字
grep '>' GRCh38.transcript.fa.tmp | head

# 去掉空格后面的字符串，保证cDNA文件中fasta序列的名字简洁，不然后续会出错
cut -f 1 -d ' ' GRCh38.transcript.fa.tmp >GRCh38.transcript.fa
head GRCh38.transcript.fa

# 定量新转录本

salmon index -t GRCh38.transcript.fa -i GRCh38.salmon

for samp in `tail -n +2 sampleFile | cut -f 1`; do \
  salmon quant --gcBias -l A -1 ${samp}_1.fq.gz -2 ${samp}_2.fq.gz \
  -i GRCh38.salmon -o ${samp}/${samp}.salmon.count -p 4 >${samp}.salmon.log 2>&1; done

for samp in `tail -n +2 sampleFile | cut -f 1`; do \
  cut -f 1,4 ${samp}/${samp}.salmon.count/quant.sf | \
  sed "1 s/TPM/${samp}/" >${samp}/${samp}.salmon.count/quant2.sf; done

paste `find . -name quant2.sf` | \
  awk 'BEGIN{OFS=FS="\t" }{line=$1; \
    for(i=2;i<=NF;i++) if(i%2==0) line=line"\t"$i; print line;}' \
  >ehbio_trans.TPM.xls
head ehbio_trans.TPM.xls

# 获取新转录本的序列

sed 's/"/\t/g' assembeCompare2Ref.annotated.gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="transcript" && $22 != "=") print $10, $22}' >assembeCompare2Ref.annotated.new.transcript.id

# 把serverData/soft目录下所有文件都拷贝到 ~/transcriptome/soft目录下
extractFastaByName.py -i GRCh38.transcript.fa -n assembeCompare2Ref.annotated.new.transcript.id >GRCh38.new_transcript.fa

grep '>' GRCh38.new_transcript.fa | head

# 编码能力预测

# 按实际修改路径
rnasamba_model=~/transcriptome/soft/full_length_weights.hdf5 

# 如果出错，运行
# pip install 'h5py==2.10.0' --force-reinstall
# iginal_keras_version = f.attrs['keras_version'].decode('utf8') AttributeError: 'str' object has no attribute 'decode'
rnasamba-classify -p GRCh38.new_transcript.translated.protein.fa GRCh38.new_transcript.codingstatus.txt GRCh38.new_transcript.fa ${rnasamba_model}

# 查看预测结果

head GRCh38.new_transcript.codingstatus.txt

# Extract non-coding RNA

grep 'noncoding' GRCh38.new_transcript.codingstatus.txt >GRCh38.new_noncoding_transcriptID.txt

extractFastaByName.py -i GRCh38.new_transcript.fa -n GRCh38.new_noncoding_transcriptID.txt >GRCh38.new_noncoding_transcript.fa

grep '>' GRCh38.new_noncoding_transcript.fa | head

grep -c '>' GRCh38.new_noncoding_transcript.fa

cut -c 1-60 GRCh38.new_noncoding_transcript.fa | head

# 查看翻译出的蛋白序列

head GRCh38.new_transcript.translated.protein.fa

# miRNA序列下载

# wget -c ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
# gunzip -c mature.fa.gz >mature.fa
# grep -A 1 '^>hsa-' genome/mature.fa | sed '/--/d' | cut -f 1 -d ' ' | sed "s/hsa-//" | sed 's/U/T/g' >genome/hsa.mature.fa


head genome/hsa.mature.fa

# miRNA-noncodingRNA靶基因预测

miranda genome/hsa.mature.fa GRCh38.new_noncoding_transcript.fa -sc 155 -en -20 | grep -B 16 '^>>' >lncRNA.miranda.predict.target

grep '>>' lncRNA.miranda.predict.target | cut -f 1,2 | sed 's/>>//' |  awk 'BEGIN{OFS=FS="\t"}{a[$1"\t"$2]+=1;}END{for(i in a) print i,a[i],"lncRNA"}' >lncRNA.miranda.predict.target.id

head lncRNA.miranda.predict.target.id

# miRNA coding RNA 3UTR 靶基因预测

miranda genome/hsa.mature.fa genome/UTR.fa -sc 155 -en -20 | grep -B 16 '^>>' >proteinCoding.miranda.predict.target

grep '>>' proteinCoding.miranda.predict.target | cut -f 1,2 | sed 's/>>//' |  awk 'BEGIN{OFS=FS="\t"}{a[$1"\t"$2]+=1;}END{for(i in a) print i,a[i],"Coding"}' >proteinCoding.miranda.predict.target.id


# miRNA sponge

awk 'BEGIN{OFS=FS="\t"; }ARGIND==1{lncTar[$1]=lncTar[$1]==""?$0:lncTar[$1]"\n"$0;}ARGIND==2{if(lncTar[$1]!="") {print $0; if(output[$1]=="") {print lncTar[$1]; output[$1]=1;}}}' lncRNA.miranda.predict.target.id proteinCoding.miranda.predict.target.id >miRNA_sponge.list.id

head miRNA_sponge.list.id

# Correlation 


awk 'BGEIN{OFS=FS="\t"}ARGIND==1{a[$2]=1;}ARGIND==2{if(FNR==1 || a[$1]==1) print $0}' miRNA_sponge.list.id ehbio_trans.TPM.xls >miRNA_sponge.list.expr

head miRNA_sponge.list.expr

pearsonCorrelation2file.py -i miRNA_sponge.list.expr -j miRNA_sponge.list.expr --positive-cor-only -o miRNA_sponge.list.cor

head miRNA_sponge.list.cor.scipy.pearson.xls

grep ENST00000246190 miRNA_sponge.list.cor.scipy.pearson.xls

selectSponge.py -i miRNA_sponge.list.id -c miRNA_sponge.list.cor.scipy.pearson.xls -o miRNA_sponge.pos_cor


head miRNA_sponge.pos_cor.table.txt

head miRNA_sponge.pos_cor.network.txt
