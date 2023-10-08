# 易微生物组(EasyMicrobiome)

易扩增子、易宏基因组等分析流程依赖常用软件、脚本文件和数据库注释文件等。

Popular software, scripts and database annotation for EasyAmplicon and EasyMetagenome

版本(Version)：EasyMicrobiome v1.20

更新时间(Update)：2023/10/13

项目主页(Project homepage): https://github.com/yongxinliu/EasyMicrobiome

## 软件安装(Install)

软件包几乎每个季度会更新生成一次，下载并添加至环境变量至可使用。

The software package will be updated and generated almost every quarter, downloaded and added to the environment variable to make it available.

多种下载软件和数据库的方法：任选其一即可
Multiple ways to download software and databases: just choose one

国内可备选微生物所下载站 http://nmdc.cn/datadownload 和百度网盘 https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315

	# 方法1. git下载，可使用wget或主页中直接下载压缩包
	git clone https://github.com/YongxinLiu/EasyMicrobiome

	# 方法2. 备用链接下载
	wget -c ftp://download.nmdc.cn/tools/soft/EasyMicrobiome.tar.gz
	tar xvzf EasyMicrobiome.tar.gz


添加linux命令可执行权限

	chmod +x EasyMicrobiome/linux/*

添加软件至环境变量，否则需要指定软件的完整路径使用

	# 临时添加环境变量
	export PATH=$PATH:`pwd`/EasyMicrobiome/linux:`pwd`/EasyMicrobiome/script"
	# 将变量写入.bashrc，永久添加环境变量
	echo "PATH=$PATH:`pwd`/EasyMicrobiome/linux:`pwd`/EasyMicrobiome/script" >> ~/.bashrc

## 使用方法

该软件为易扩增子、易宏基因组的依赖包。详细使用见各项目主页：

- 易扩增子分析流程：https://github.com/YongxinLiu/EasyAmplicon

- 易宏基因组分析流程：https://github.com/YongxinLiu/EasyMetagenome

流程的绘图部分，依赖的R包较多，推荐在Windows系统是使用(安装R包更方便)，同时提供了4百个包的合集下载，节省安装时间

- R语言4.3环境和R包：R语言主页 http://www.r-project.org ，Windows版包合集 ftp://download.nmdc.cn/tools/win/4.3.zip

## 软件清单

*注：名称的链接对应软件的主页，大部分已经整合入本项目。对于较大的文件，标题后提供下载链接，使用时需自行下载。

- linux：Linux系统下分析软件
    - [microbiome_helper](https://github.com/LangilleLab/microbiome_helper)：微生物组分析输助脚本，如metaphlan2结果转换STAMP格式(metaphlan_to_stamp.pl)，picurst结果功能组成绘图(plot_metagenome_contributions.R)
    - Miniconda3-latest-Linux-x86_64.sh：软件管理器 https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    - qiime2-2023.2.tar.gz：QIIME2安装包，解压至conda的envs目录可用 ftp://download.nmdc.cn/tools/conda/qiime2-2023.2.tar.gz
    - qiime2-2023.2-py38-linux-conda.yml：QIIME2软件安装清单，使用conda在线安装
    - [sparcc](https://github.com/TankMermaid/sparcc)：sparcc网络分析python脚本
    - [usearch](http://www.drive5.com/usearch/)：扩增子分析流程
    - [vsearch](https://github.com/torognes/vsearch)：扩增子分析流程(免费64位版usearch)
- mac：Mac系统下分析软件
    - csvtk：表格分析工具
    - iqtree：进化树构建
    - qiime2-2023.2-py38-osx-conda.yml：QIIME2软件安装清单，使用conda在线安装
    - R-4.2.3.pkg：R语言安装包
    - RStudio-2023.03.0-386.dmg：RStudio安装包
    - rush：并行管理工具
    - seqkit：序列处理工具
    - taxonkit：NCBI分类处理工具
    - usearch：扩增子分析流程
    - vsearch：扩增子分析流程(免费64位版usearch)
- win：Windows系统下分析软件
    - [Git-2.40.0-64-bit.exe](http://gitforwindows.org/)：提供Git bash环境，自行下载安装，教程见：[Windows轻松实现linux shell环境：gitforwindows](https://mp.weixin.qq.com/s/KtM4c4o4iLfD4ZkEnMi1pg)
    - [R-4.2.3-win.exe](https://www.r-project.org/ )：R语言安装包，下载最新版：Downad CRAN - China Tsinghua - Download R for Windows(Mac) —— base —— Download R 4.2.0
    - [RStudio-2023.03.0-386.exe](https://www.rstudio.com/products/rstudio/download/#download)：RStudio安装包，提供分析运行界面。
    - [4.2.zip](ftp://download.nmdc.cn/tools/win/4.2.zip)：R语言常用400+包合集，解压至R包安装位置即可用。
    - [usearch.exe](http://www.drive5.com/usearch/)：扩增子分析流程
    - [vsearch.exe](https://github.com/torognes/vsearch)：扩增子分析流程(免费64位版usearch)
    - [STAMP2.1.3](http://kiwi.cs.dal.ca/Software/STAMP)：微生物组图形界面差异分析工具
    - Adobe_Illustrator_CC_2018_v22.1.0.314_x64_zh_CN_Portable.7z：图片拼图、模式图绘制工具，使用试用版或自行购买
    - [Cytoscape_3_8_2_windows_64bit.exe](http://www.cytoscape.org)：网络分析安装包
    - [csvtk.exe](https://github.com/shenwei356/csvtk)：表格分析工具
    - [seqkit.exe](https://github.com/shenwei356/seqkit)：序列处理工具
    - [taxonkit.exe](https://github.com/shenwei356/taxonkit )：NCBI分类处理工具
    - [rush.exe](https://github.com/shenwei356/rush)：并行管理工具
    - [epp510_1828_64bit.exe](https://www.editplus.com/)：文本编辑器
    - [Xshell](https://www.netsarang.com/zh/free-for-home-school)：远程访问服务器终端，需要申请免费版下载链接；备选[PuTTY](http://www.putty.be/)
    - [FileZilla](https://filezilla-project.org/download.php?type=client)：远程访问服务器文件上传下载，备选[WinSCP](https://winscp.net/eng/download.php)
    - [gephi-0.9.2-windows.exe](https://gephi.org/)：网络图绘制工具
    - iqtree.exe：进化树构建
    - libiomp5md.dll：动态库，iqtree运行中提示缺少时，可添加至软件所在目录
    - jdk-11.0.7_windows-x64_bin.exe：Java运行环境
    - muscle.exe：多序列比对工具
    - npp.7.8.9.Installer.x64.exe：文本编辑器NotePad++安装包
    - rtools40-x86_64.exe：R源码安装时的编绎工具
    - wget.exe：命令行下载工具

## 数据库

- gg：GreenGenes细菌16S数据库
    - [gg_13_8_otus.tar.gz](ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
)：13年8月更新OTU数据库，用于usearch有参定量和PICRUSt/BugBase功能预测、QIIME 2制作分类器。[国内备份链接](ftp://download.nmdc.cn/tools/gg/gg_13_8_otus.tar.gz)
    - [16S_13_5_precalculated.tab.gz](ftp://download.nmdc.cn/tools/picrust/16S_13_5_precalculated.tab.gz)：picrust的GreenGenes 16S拷贝数
	- [ko_13_5_precalculated.tab.gz](ftp://download.nmdc.cn/tools/picrust/ko_13_5_precalculated.tab.gz)：picrust的GreenGenes 16S对应的KO数量信息
- kegg：KEGG数据库描述信息整理
    - [ko00001.keg](https://www.kegg.jp)：KEGG层级注释体系，主页 —— KEGG BRITE —— KEGG Orthology (KO) —— Download htext，下载保存为ko00001.tsv
    - ko00001.tsv：转换jason格式为制表符分隔的KO对应描述、(三级)通路、二级通路和一级通路信息
    - KO1-4.txt：KO对应的3级注释，包括(三级)通路、二级通路和一级通路信息，用于KO表的分类汇总
    - KO_description.txt：KO编号对应的功能描述
    - KO_path.list：KO与通路(Pathway)的对应关系，存在某个KO存在于多个通路(1对多)
- usearch：usearch/vsearch物种分类sintax命令使用数据库
    - [rdp_16s_v18.fa.gz](http://www.drive5.com/usearch/manual/sintax_downloads.html)：16S的RDP16数据库，usearch作者整理，更多16S、ITS和18S数据库见 http://www.drive5.com/usearch/manual/sintax_downloads.html
    - rdp_16s_v18.fa.gz：16S的RDP18数据库，2021年基于RDP数据库整理
    - utax_reference_dataset_all_04.02.2020.fasta.gz：ITS注释数据库，可从UNITE下载
- eggnog: eggnog结果的注释文件补充
    - COG.anno：COG的第一、二级注释

## 脚本 

- 使用说明：分析常用脚本类型
    - .R文件为R脚本，使用Rscript命令执行；
    - .sh为Shell脚本，使用/bin/bash命令执行；
    - .pl为Perl脚本，使用perl命令执行；
    - .py为Python脚本，使用python执行，注意还分为python2和python3两种

- script：微生物组数据分析
    - BugBase：16S扩增子表型预测R脚本和数据库
    - FAPROTAX_1.2.4：16S扩增子元素循环预测Python脚本和数据库
    - table2itol：iTOL进化树注释文件制作R脚本
    - alpha_barplot.R：Alpha多样性指数柱状图+标准差图绘制
    - alpha_boxplot.R：Alpha多样性指数箱线图+统计绘制
    - alpha_rare_curve.R：usearch计算稀释曲线可视化
    - beta_cpcoa.R：基于距离矩阵开展限制性PCoA分析及可视化散点图+分组着色+置信椭圆，要求至少3个分组
    - beta_pcoa.R：基于距离矩阵的主坐标PCoA分析及可视化散点图+分组着色+置信椭圆+组间两两统计
    - BetaDiv.R：更多Beta多样性分析，如PCA、PCoA、NMDS、LDA、CCA、RDA等
    - compare.R：两组比较，支持t.test、wilcox、edgeR三种方法
    - compare_heatmap.R/sh：基于两组比较结果绘制热图
    - compare_manhattan.sh：基于两组比较结果绘制曼哈顿图
    - compare_volcano.R：基于两组比较结果绘制火山图
    - faprotax_report_sum.pl：FARPROTAX分析结果报告整理
    - filter_feature_table.R：按频率过滤OTU表
    - format_dbcan2list.pl：dbcan数据库注释结果整理
    - format2lefse.R：OTU表和物种注释生成LEfSe输入文件
    - format2stamp.R：OTU表和物种注释生成STAMP输入文件
    - kegg_ko00001_htext2tsv.pl：KEGG注释结果整理
    - kraken2alpha.R：Kraken2结果整理、抽平和alpha多样性指数计算
    - mat_gene2ko.R：按类型折叠表格
    - metaphlan_boxplot.R：metaphalan2结果可视化为箱线图
    - metaphlan_hclust_heatmap.R：metaphalan2结果可视化为聚类热图
    - metaphlan_to_stamp.pl：metaphalan2结果转换为STAMP格式
    - otu_mean.R：OTU表统计分组均值(总体均值)、分组求合
    - otutab_filter_nonBac.R：16S的OTU表按sintax注释结果选择细菌、古菌且过滤叶绿体和线粒体
    - otutab_filter_nonFungi.R：ITS的OTU表选择真菌
    - otutab_freq2count.R：转换频率为伪整数，用于要求整型输入的分析，如多样性、edgeR差异分析等
    - otutab_rare.R：OTU表抽平
    - plot_metagenome_contributions.R：PICRUSt结果物种的功能组成绘制
    - sp_pheatmap.sh：绘制热图
    - sp_vennDiagram.sh：绘制维恩图
    - summarizeAbundance.py：按类型折叠大表，如基因按KEGG的KO合并
    - tax_circlize.R：物种组成圈图
    - tax_maptree.R：物种组成气泡图
    - tax_stackplot.R：物种组成堆叠柱状图

使用此脚本，请引用下文：

If used this script, please cited:

**Yong-Xin Liu**, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin Zhou, Liang Chen, Xubo Qian, Jiao Xi, Hongye Lu, Huiluo Cao, Xiaoya Ma, Bian Bian, Pengfan Zhang, Jiqiu Wu, Ren-You Gan, Baolei Jia, Linyang Sun, Zhicheng Ju, Yunyun Gao, **Tao Wen**, **Tong Chen**. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. **iMeta** 2: e83. https://doi.org/10.1002/imt2.83

Copyright 2016-2023 Yong-Xin Liu <liuyongxin@caas.cn>, Tao Wen <taowen@njau.edu.cn>, Tong Chen <chent@nrc.ac.cn>


