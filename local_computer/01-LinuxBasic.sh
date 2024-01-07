# Last login: Mon Jun  5 16:56:56 2017 from 239.241.208.209
# Welcome to aliyun Elastic Compute Service!
# ct@ehbio:~$


# 若提示命令找不到，运行下面语句安装tree
# 需要有根用户权限
# yum install tree.x86_64
tree -d -L 2 /

df -h

# Filesystem            Size  Used Avail Use% Mounted on
# /dev/sda2             193G   61G  122G  34% /
# tmpfs                 127G  344K  127G   1% /dev/shm
# /dev/sda1             190M   77M  103M  43% /boot
# /dev/mapper/a          37T   12T   25T  32% /ehbio1
# /dev/mapper/b          37T   28T  8.8T  76% /ehbio2
# /dev/mapper/c          37T   15T   23T  40% /ehbio3


# serverInfo.sh是我写的一个脚本，这个脚本怎么实现的会是一个考核题目。
serverInfo.sh

ls

is
#-bash: is: 未找到命令
# 大小写敏感

lS
#-bash: lS: 未找到命令


mkdir data
ls
data

mkdir data

#mkdir: 无法创建目录"data" : 文件已存在

# -p: no error if existing, make parent directories as needed
mkdir -p data


mkdir data
cat <<END
a
bc
END


cat <<END >data/ehbio.fa
>SOX2
ACGTCGGCGGAGGGTGGSCGGGGGGGGAGAGGT
ACGATGAGGAGTAGGAGAGAGGAGG
>OCT4
ACGTAGGATGGAGGAGAGGGAGGGGGGAGGAGAGGAA
AGAGTAGAGAGA
>NANOG
ACGATGCGATGCAGCGTTTTTTTTTGGTTGGATCT
CAGGTAGGAGCGAGGAGGCAGCGGCGGATGCAGGCA
ACGGTAGCGAGTC
>mYC HAHA
ACGGAGCGAGCTAGTGCAGCGAGGAGCTGAGTCGAGC
CAGGACAGGAGCTA
end
END

## 注意命令和参数之间的空格
ls-l

#-bash: ls-l: 未找到命令

ls -l

#总用量 4
## d: dir; 表示data是个目录
## rwx：表示目录的权限，暂时忽略，或自己在线搜索
#drwxrwxr-x 2 ct ct 4096 6月   8 14:52 data

ls -l data
#总用量 4
## 开头的`-`表示ehbio.fa是个文件
#-rw-rw-r-- 1 ct ct 284 6月   8 14:48 ehbio.fa


# 这个应该是最常见的错误之一，程序不可能知道你的输入文件在什么地方，
# 需要人为指定。
# 如果未指定路径，表示当前目录
cat ehbio.fa
cat: ehbio.fa: 没有那个文件或目录

cat data/ehbio.fa


# 记住输出
pwd

#/home/ct

cd data

# 注意输出变化
pwd
/home/ct/data

head -n 6 ehbio.fa


less ehbio.fa
# q: 退出
# 上下箭头、空格翻页


man ls


# gzip -c 把压缩的文件输出到标准输出 (一般是屏幕)
# '>' 输出重定向，输出写入文件

gzip -c ehbio.fa >ehbio.fa.gz

# 多了一个.gz文件
ls
#ehbio3.fa  ehbio4.fa  ehbio5.fa  ehbio.fa  ehbio.fa.gz  second.fa

# 解压缩
gunzip ehbio.fa.gz
gzip: ehbio.fa already exists; do you wish to overwrite (y or n)? y

ls
#ehbio3.fa  ehbio4.fa  ehbio5.fa  ehbio.fa  second.fa


# 输出文件有14行
wc -l ehbio.fa

#14 ehbio.fa


cut -f 1 -d ' ' ehbio.fa | tail -n 4

# >mYC
# ACGGAGCGAGCTAGTGCAGCGAGGAGCTGAGTCGAGC
# CAGGACAGGAGCTA
# end
man cut

# 直接跳到上面运行的cut命令，再执行一次
!cut
cut -f 1 -d ' ' ehbio.fa | tail -n 4
# >mYC
# ACGGAGCGAGCTAGTGCAGCGAGGAGCTGAGTCGAGC
# CAGGACAGGAGCTA
# end


# 写完下面的命令，突然不想运行了，又不想一个个删掉
cut -f 1 -d ' ' ehbio.fa | tail -n 4

# 按ctrl+a, 回到行首，再输入`#`号，回车，命令即被注释掉。
#cut -f 1 -d ' ' ehbio.fa | tail -n 4



# 管道符的使用
# 第一个命令的输出作为第二个的输入
# 前面的例子中也有使用
# tr: 是用于替换字符的，把空格替换为换行，文字就从一行变为了一列
echo "1 2 3" | tr ' ' '\n'

# 1
# 2
# 3



ls -l /home
#-rwxr-xr-x 1 ct ct 26 12月  7 2016  ct

# 让自己的家目录只自己可见
chmod go-rwx /home/ct
ls -l /home
#-rwx------ 1 ct ct 26 12月  7 2016  ct

# 同组人增加可读和可执行属性
chmod g+rx /home/ct
ls -l /home
-rwxr-x--- 1 ct ct 26 12月  7 2016  ct

# 所有人增加可读和可执行属性
chmod a+rx /home/ct
ls -l /home
-rwxr-xr-x 1 ct ct 26 12月  7 2016  ct


ls -l "`which cd`"
#rwx: 文件所有者可读、可写、可执行
#r-x: 文件所有者所在组其它成员可读、可执行，不可修改
#r-x: 其它人可读、可执行，不可修改
#-rwxr-xr-x 1 root root 26 12月  7 2016 /usr/bin/cd

ls -l "`which mkdir`"
#-rwxr-xr-x. 1 root root 79768 11月  6 2016 /usr/bin/mkdir

ls -l "`which python`"
#l: 代表软连接
#软连接自身是所有人可读可写，但具体的权限依赖于其链接的文件
#lrwxrwxrwx. 1 root root 7 3月  22 15:04 /usr/bin/python -> python2


# 新建个文件
cat <<END >run.sh
echo " I am a script created by ehbio."
END

# 查看其权限值
ls -l run.sh
# -rw-rw-r-- 1 ct ct 39 6月  14 23:12 run.sh

# 更改权限值
chmod 755 run.sh

# 查看其权限值
# 注意多了3个x
ls -l run.sh
# -rwxr-xr-x 1 ct ct 39 6月  14 23:12 run.sh

# 去除其它用户的可执行权限
chmod o-x run.sh

# 注意看少了个x
ls -l run.sh
# -rwxr-xr-- 1 ct ct 39 6月  14 23:12 run.sh

# 去除同组的可执行权限
chmod g-x run.sh

# 注意看又少了个x
ls -l run.sh
# -rwxr--r-- 1 ct ct 39 6月  14 23:12 run.sh

# 去除所有人的可执行权限
chmod a-x run.sh
# ls -l run.sh
# -rw-r--r-- 1 ct ct 39 6月  14 23:12 run.sh

# 给所有人增加可执行权限
chmod a+x run.sh
# ls -l run.sh
# -rwxr-xr-x 1 ct ct 39 6月  14 23:12 run.sh


run.sh
#-bash: run.sh: 未找到命令


echo $PATH
#/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin


# 加到环境变量的路径必须是全路径，全路径指以/开头或已~开头的路径
# 注意第一个PATH不含$, 第二个PATH有$符号
# 我们后面会讲什么时候用$, 什么时候不用$
export PATH=$PATH:/home/ct
echo $PATH
/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/ct


run.sh
I am a script created by ehbio.


# 这是我的~/.bash_profile中的内容，主要是最好一行。可以连续的加入多个路径。
if [ -f ~/.bashrc ]; then
. ~/.bashrc
fi

if [ -f ~/.bash_aliases ]; then
. ~/.bash_aliases
fi

export PATH=$PATH:/home/ct:/home/bin:/home/soft/bowtie2/bin


# 注意$PATH的顺序

export PATH=/home/ct/anaconda/bin:$PATH


