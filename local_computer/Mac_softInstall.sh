mkdir -p ${HOME}/bin

cat <<END >~/.bashrc
#echo "export PATH=/usr/local/opt/coreutils/libexec/gnubin:\$PATH:${HOME}/bin"
alias awk=gawk
alias sed=gsed
END

source ~/.bashrc

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
#brew install bash
brew install coreutils

brew install gawk
brew install gnu-sed

# 程序都在，且都是绿色为正确
ls -l ${HOME}/bin/*

cd public/mac
unzip R4.1_mac.packages.zip
cp -rf Library/* /Library/Frameworks/R.framework/Versions/4.1/Resources/library/

# 下载vsearch

# wget https://github.com/torognes/vsearch/releases/download/v2.14.1/vsearch-2.14.1-macos-x86_64.tar.gz
# tar xvzf vsearch-2.14.1-macos-x86_64.tar.gz
# mv vsearch-2.14.1-macos-x86_64/bin/vsearch ~/bin

# 安装Muscle

wget https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86darwin64.tar.gz
tar vxzf muscle3.8.31_i86darwin64.tar.gz
mv muscle3.8.31_i86darwin64 ~/bin/muscle

# IQtree

wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-MacOSX.zip
unzip iqtree-1.6.12-MacOSX.zip
mv iqtree-1.6.12-MacOSX/bin/iqtree ~/bin





# cp public/mac/usearch ${HOME}/bin
# cp public/mac/vsearch ${HOME}/bin

# cd 02Software/mac/
# unzip iqtree-1.6.10-MacOSX.zip

# mv iqtree-1.6.10-MacOSX/bin/iqtree ${HOME}/bin

chmod 755 ${HOME}/bin/*

# 程序都在，且都是绿色为正确
ls -l ${HOME}/bin/*




