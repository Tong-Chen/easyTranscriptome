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


## 转录本拼装

使用STAR比对的结果拼装时，一定要加比对参数`--outSAMattrIHstart 0 --outSAMstrandField intronMotif`，不然出来的都是单外显子转录本。

#单样本拼装


# # -G: 指定reference GTF
# -p 1: 使用一个线程，多核处理器可调大
# -f 0.01 : 允许的最小isoform比例，默认0.01

cd ${wd}
stringtie trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam \
  -G genome/GRCh38.gtf -l trt_N061011 -o trt_N061011/trt_N061011.stringtie_first.gtf \
  -f 0.01 -p 2
  
head trt_N061011/trt_N061011.stringtie_first.gtf
  
# 多样本循环拼装

cat sampleFile

cd ${wd}
for i in `tail -n +2 sampleFile | cut -f 1`; do 
	stringtie ${i}/${i}.Aligned.sortedByCoord.out.bam -G genome/GRCh38.gtf -l ${i} -o ${i}/${i}.stringtie_first.gtf -f 0.01 -p 1
done &



# 转录本合并
cd ${wd}
# 获取所有拼装好的gtf文件
find . -name *.stringtie_first.gtf >mergeList.txt

# -G: 指定reference GTF
# -l: 输出结果中转录本的名字前缀
# -o: 输出文件
# mergeList.txt：单个样品GTF列表，每个一行

stringtie --merge -G genome/GRCh38.gtf -l ehbio_trans -o ehbio_trans.gtf mergeList.txt


#新拼装转录本与原基因组注释转录本比较 (可用来筛选新转录本)

# -R: 只考虑与拼装的转录本有重叠的注释转录本
# -r: reference gtf
# -o: 输出前缀

cd ${wd}
gffcompare -R -r genome/GRCh38.gtf -o assembeCompare2Ref ehbio_trans.gtf
# 会输出一个assembeCompare2Ref.annotated.gtf，用于后续的定量

head assembeCompare2Ref.annotated.gtf

统计不同种类转录本的数目

cut -f 3 assembeCompare2Ref.ehbio_trans.gtf.tmap | tail -n +2 | sort | uniq -c

# 转录本定量

## 基于HTSEQ

for i in `tail -n +2 sampleFile | cut -f 1`; do 
	htseq-count -f bam -r pos -a 10 -t exon -s no -i gene_id -m union ${i}/${i}.Aligned.sortedByCoord.out.bam assembeCompare2Ref.annotated.gtf >${i}/${i}.readsCount
	grep -v '^__' ${i}/${i}.readsCount | sed "1 iGene\t${i}"  >${i}/${i}.readsCount2
done &

## 基于Salmon



## 自己完成

# 合并表达文件和差异分析同使用基因组注释文件计算差异基因一致。

cd ${wd}
paste `find . -name *.readsCount2` | \
  awk 'BEGIN{OFS=FS="\t" }{line=$1; \
    for(i=2;i<=NF;i++) if(i%2==0) line=line"\t"$i; print line;}' \
  >ehbio_trans.Count_matrix2.xls
head ehbio_trans.Count_matrix2.xls

# 差异剪接分析 {#alternative_splicing}


# -b1, -b2参数的值是一个文件，文件内是一个样品多个重复的bam文件
# 所有bam文件写在一行，用逗号隔开
# 下面例子中考虑到每行能展示的宽度有限，所以文件内容做了换行处理，实际是一行
cd ${wd}

awk 'BEGIN{OFS=FS="\t"}{if(FNR>1) a[$2]=a[$2]==""?"./"$1"/"$1".Aligned.sortedByCoord.out.bam":a[$2]",./"$1"/"$1".Aligned.sortedByCoord.out.bam"}END{for(i in a) {file_name=i".bam.txt"; print a[i] >file_name;}}' sampleFile

# cat <<END >trt.bam.txt
# ./trt_N052611/trt_N052611.Aligned.sortedByCoord.out.bam,./trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam,./trt_N080611/trt_N080611.Aligned.sortedByCoord.out.bam,./trt_N61311/trt_N61311.Aligned.sortedByCoord.out.bam
# END

# cat <<END >untrt.bam.txt
# ./untrt_N052611/untrt_N052611.Aligned.sortedByCoord.out.bam,./untrt_N061011/untrt_N061011.Aligned.sortedByCoord.out.bam,./untrt_N080611/untrt_N080611.Aligned.sortedByCoord.out.bam,./untrt_N61311/untrt_N61311.Aligned.sortedByCoord.out.bam
# END

# --gtf: 指定基因注释文件，可以是自己下载的注释，也可以是stringtie拼装完merge后的注释
# --od: 输出目录 (output dir)
# -t: paired-end or single-end
# --libtype: 链特异性类型
# -nthread: 多线程计算; --tstat: 统计模型多线程计算
# --cstat: 差异剪接阈值，默认0.0001, 表示有0.01%的差异，取值在0-1之间

rmats.py --b1 trt.bam.txt --b2 untrt.bam.txt --gtf genome/GRCh38.gtf --od trt_untrt -t paired --libType fr-unstranded --readLength 63 --nthread 2 --tstat 2 --cstat 0.0001 --tmp tmp


#筛选SE差异显著的剪接位点

awk 'FNR==1 || $20<0.2' trt_untrt/SE.MATS.JC.txt >trt_untrt/SE.MATS.JC.sig.txt

head trt_untrt/SE.MATS.JC.sig.txt

#筛选MXE差异显著的剪接位点

awk 'FNR==1 || $20<0.2' trt_untrt/MXE.MATS.JC.txt >trt_untrt/MXE.MATS.JC.sig.txt

head trt_untrt/MXE.MATS.JC.sig.txt

#  sashimiplot

# --b1, --b2 同rMATS，只是直接跟文件内容而不是文件名
# -t: 想要绘制的类型
# -e: rMATS的输出，一般选择差异显著的绘制
# --exon_s,  --intron_s: 缩放外显子或内含子，默认无缩放；一般用于内含子太大时，把内含子
# 相对缩小一点，图会更好看一些
# --l1, --l2，对应于--b1, --b2的样品的组名
# -o 指定输出目录

rmats2sashimiplot --b1 ./trt_N052611/trt_N052611.Aligned.sortedByCoord.out.bam,./trt_N061011/trt_N061011.Aligned.sortedByCoord.out.bam,./trt_N080611/trt_N080611.Aligned.sortedByCoord.out.bam,./trt_N61311/trt_N61311.Aligned.sortedByCoord.out.bam \
  --b2 ./untrt_N052611/untrt_N052611.Aligned.sortedByCoord.out.bam,./untrt_N061011/untrt_N061011.Aligned.sortedByCoord.out.bam,./untrt_N080611/untrt_N080611.Aligned.sortedByCoord.out.bam,./untrt_N61311/untrt_N61311.Aligned.sortedByCoord.out.bam \
  -t SE -e trt_untrt/SE.MATS.JC.sig.txt --exon_s 1 --intron_s 5 --l1 trt --l2 untrt -o trt_untrt/SE.MATS.JS.sig.sashimiplot


