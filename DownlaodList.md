[TOC]

## 微生物组软件包

### 易微生物组(EasyMicrobiome)脚本库

主页：https://github.com/YongxinLiu/EasyMicrobiome

简介：易扩增子、易宏基因组等分析流程依赖常用软件、脚本文件和数据库注释文件等。

版本 | 大小 | 更新时间 | 下载链接 | 描述
-|-|-|-|-
EasyMicrobiome v1.18 | 174M | 2023-04-04 | [zip](ftp://download.nmdc.cn/tools/soft/EasyMicrobiome.tar.gz) | 分析脚本、命令行程序、数据库注释整理表格等
Windows系统常用软件 | 2G | 2023-04-04 | [目录链接](ftp://download.nmdc.cn/tools/win/win.zip) | 常用软件Git、R、RStudio、Rtools、STAMP
微生物组常用软件补充目录 | 2G | 2023-03-22 | [目录链接](ftp://download.nmdc.cn/tools/soft/) | 常用软件、数据库等


### R语言安装包

包有585多个常R包的安装版。解决安装中无法下载、编绎失败、缺少依赖关系等问题。下载内容解压后替换R包安装目录。R包位置查看：在RStudio中右下子窗口 Packages菜点，点击Install，窗口中Install to Library即为安装位置，如我的Win11系统为 C:/User/yongxinliu/AppData/Local/R/win-library/4.2

版本 | 大小 | 更新时间 | 下载链接
-|-|-|-
Windows10版4.2.x | 942MB | 2023-03-30 | [zip](ftp://download.nmdc.cn/tools/win/4.2.zip)
MacOS版4.2.x | 606MB | 2023-03-30 | [zip](ftp://download.nmdc.cn/tools/mac/R4.2_mac_libraryX86_64.zip)


### 易扩增子(EasyAmplicon)分析流程

主页：https://github.com/YongxinLiu/EasyAmplicon

简介：跨平台(Win/Mac/Linux)的扩增子分析流程，在RStudio中完成可重复分析，去噪/分类注释比QIIME 2快上百倍，提供从原始数据到特征表，以及数十种出版级别的可视化方案，详见iMeta文章[英文版](https://onlinelibrary.wiley.com/doi/10.1002/imt2.83)/[中文版](https://mp.weixin.qq.com/s/oa7QDDNKLU1TmQFPFhcPGQ)。

版本 | 大小 | 更新时间 | 下载链接 | 描述
-|-|-|-|-
EasyAmplicon v1.18 | 180M | 2023-03-11 | [zip](ftp://download.nmdc.cn/tools/soft/EasyAmplicon.tar.gz) | 流程的脚本、依赖命令行程序、RDP/UNITE数据库和说明文档

### 易宏基因组(EasyMetagenome)分析流程

主页：https://github.com/YongxinLiu/EasyMetagenome

简介：易宏基因组分析流程代码，以及示例数据分析结果。需提交安装易微生物组才可以运行。

版本 | 大小 | 更新时间 | 下载链接 | 描述
-|-|-|-|-
EasyMetagenome v1.18 | 121M | 2023-04-04 | [zip](ftp://download.nmdc.cn/tools/soft/EasyMetagenome.tar.gz) | 分析流程脚本、结果示例等

### Conda软件压缩包

Conda虽然极大方便了软件的安装，但使用中经常出现环境冲突安装失败的情况。这里提供安装好的Conda环境打包，可以下载解压实现快速安装，并极大提高成功率。

版本 | 大小 | 更新时间 | 下载链接 | 描述
-|-|-|-|-
QIIME2 2023.2 | 1.5G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/qiime2-2023.2.tar.gz) | 扩增子分析流程
picrust | 232M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/picrust.tar.gz) | 功能注释
picrust2 | 490M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/picrust2.tar.gz) | 功能注释新版，要求16G内存起
[KneadData v0.12.0](https://bioconda.github.io/recipes/kneaddata/README.html) | 643M | 2023-04-01 | [tar.gz](ftp://download.nmdc.cn/tools/conda/kneaddata.tar.gz) | 宏基因组序列质量控制和去宿主流程
HUMAnN2+MetaPhlAn2 | 400M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/humann2.tar.gz) | HUMAnN 2物种和功能注释流程，包括humann2 v2.8.1、metaphlan2 2.7.5、graphlan 1.1.3、export2graphlan 0.22等
LEfSe 1.1.2 | 441M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/lefse.tar.gz) | 生物标志物鉴定，结果可绘制为柱状图、物种树等
Kraken2.1.2 | 527M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/kraken2.tar.gz) | 打包的kraken2流程，包括kraken2、bracken、krakentools、krona、r-optparse等
megahit/spades组装流程 | 768M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/megahit.tar.gz) | 宏基因组组装流程，包括megahit v1.2.9、metaSPAdes v3.15.4、MetaQUAST v5.0.2、CD-HIT v4.8.1、EMBOSS v6.6、salmon v1.8等软件
eggnog-mapper 2.1.10 | 227M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/eggnog.tar.gz) | 免费全面的功能注释软件，可实现KO/COG/CAZy/GO等10余种注释信息，需要配合数据库使用，KEGG注释的替换工具
RGI 5.2.1 | 329M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/rgi.tar.gz) | 抗生素抗性/耐药基因注释软件
MetaWRAP 1.3 | 1.7G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/metawrap.tar.gz) | 分箱流程，包括组装、提纯、定量
dRep v3.2.2 | 246M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/drep.tar.gz) | 细菌基因组、宏基因组组装基因组去冗余
GTDB-tk v2.2.6 | 154M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/gtdbtk.tar.gz) | 细菌、古菌基因组物种注释

## 扩增子数据库

### EzBioCloud

主页：https://www.ezbiocloud.net/

简介：EzBioCloud是综合的细菌16S鉴定数据库，所有16S序列经人工校正，几乎全部为完整27F-1492R全长16S序列，而且全面覆盖NCBI、JGI的16S和细菌基因组，以及PacBio测序的16S全长序列。数据库每季度更新，近10年来被引用过万次。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
USEARCH | 19M | 2021-07-07 | [fasta.gz](ftp://download.nmdc.cn/tools/amplicon/EzBioCloud/ezbiocloud_usearch.fasta.gz) | usearch/vsearch格式物种注释数据库
MOTHUR | 19M | 2021-07-07 | [fasta.gz](ftp://download.nmdc.cn/tools/amplicon/EzBioCloud/EzBioCloud_16S_database_for_MOTHUR.zip) | MOTHUR 格式物种注释数据库
QIIME | 19M | 2021-07-07 | [fasta.gz](ftp://download.nmdc.cn/tools/amplicon/EzBioCloud/EzBioCloud_16S_database_for_QIIME.zip) | usearch/vsearch格式物种注释数据库

### GreenGenes

主页：http://ftp.microbio.me/greengenes_release/current/

简介：GreenGenes是使用最广泛的16S rRNA基因注释数据库，QIIME/QIIME 2的默认数据库，也是使用PICRUSt/BugBase等功能预测的参考数据库。2022年10月更新。QIIME 2版下载链接https://docs.qiime2.org/2023.2/data-resources/

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
2022.10.taxonomy.asv.nwk.qza | 179M | 2023-03-30 |  [qza](ftp://download.nmdc.cn/tools/amplicon/GreenGenes/2022.10.backbone.full-length.nb.qza) | QIIME2使用的物种注释文件，在QIIME 2 2023.3版本上测试通过
gg_13_8_otus | 320M | 2013-08-30 |  [tar.gz](ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus.tar.gz) | 合集，包括97/99%等多种相似度的序列、物种注释和多序列对齐文件
gg_13_8_otus_97 | 29M | 2013-08-30 | [fasta.gz](ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus_97.tar.gz) | 97%聚类序列，PICRUSt和BugBase分析准备输入文件时参考数据库

### RDP

主页：http://rdp.cme.msu.edu/

简介：USEARCH推荐使用的物种注释数据库，其精选的训练集准确率高，体积小巧速度快。RDP网站还有在线扩增子分析流程rdpipeline和引用过万的分类器(RDP Classifier)。USEARCH提供了各版本的物种注释数据库 https://www.drive5.com/usearch/manual/sintax_downloads.html 。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
RDP 16s v18 | 12.0M | 2021-09-15 | [zip](ftp://download.nmdc.cn/tools/amplicon/RDP/RDPClassifier_16S_trainsetNo18_rawtrainingdata.zip) | RDP物种注释数据库
QIIME rdp_16s_v18 | 6.0M | 2021-09-15 | [zip](ftp://download.nmdc.cn/tools/amplicon/RDP/RDPClassifier_16S_trainsetNo18_QiimeFormat.zip) | QIIME物种注释数据库
usearch rdp_16s_v18 | 6.8M | 2021-09-15  | [fa.gz](ftp://download.nmdc.cn/tools/amplicon/usearch/rdp_16s_v18.fa.gz) | usearch/vsearch物种注释数据库，推荐属水平
usearch rdp_16s_v18_sp | 6.0M | 2021-09-15 | [fa.gz](ftp://download.nmdc.cn/tools/amplicon/usearch/rdp_16s_v18_sp.fa.gz) | usearch/vsearch物种注释数据库，自制种水平

### SILVA

主页：https://www.arb-silva.de/

简介：最新最全的16S/18S核糖体数据库，最近更新于2020年 https://www.arb-silva.de/no_cache/download/archive/current/Exports/，QIIME和USEARCH都可以相应格式的数据库可用。但数据库较大，使用会消耗更多的计算资源和时间。usearch版数据库下载：https://www.drive5.com/usearch/manual/sintax_downloads.html；QIIME 2版数据下载

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
silva-138-99 | 531M | 2023-04-02 | [qza](ftp://download.nmdc.cn/tools/amplicon/silva/silva-138-99-nb-classifier.qza) | QIIME 2 2023.2物种注释数据库
silva_16s_v123 | 438M | 2020-09-23 | [fa.gz](ftp://download.nmdc.cn/tools/amplicon/silva/silva_16s_v123.fa.gz) | usearch/vsearch物种注释数据库

### UNITE

主页：https://unite.ut.ee/

简介：最好的真菌/真核微生物转录间格区(ITS)数据库。支持主流工具，如QIIME/USEARCH/Mothur等。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
UNITE 9.0 usearch格式 | 39M | 2022-11-29 | [fasta.gz](ftp://download.nmdc.cn/tools/amplicon/usearch/utax_reference_dataset_all_29.11.2022.fasta.gz) | usearch/vsearch真核生物包括真菌物种注释数据库



## 宏基因组软件和数据库

### KneadData质控去宿主

主页：http://huttenhower.sph.harvard.edu/kneaddata

简介：用于宏基因组质控、去宿主、宏转录组去核糖体等功能。下面提供常用数据库的Bowtie2索引。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
kneaddata 0.12.0 | 673M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/kneaddata.tar.gz) | 宏基因组质控、去宿主软件流程，包括kneaddata 0.12.0、fstqc v0.12.1、multiqc 1.13等
人基因组hg37 | 3.7G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/kneaddata/human_genome/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz) | 鉴定人类样本宿主含量或人源污染

### MetaPhlAn+HUMAnN物种和功能注释

- HUMAnN2+MetaPhlAn2数据库

主页：http://www.huttenhower.org/humann2

简介：宏基因组数据有参快速物种和功能通路定量软件。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
HUMAnN2+MetaPhlAn2 | 400M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/humann2.tar.gz) | HUMAnN 2物种和功能注释流程，包括humann2 v2.8.1、metaphlan2 2.7.5、graphlan 1.1.3、export2graphlan 0.22等
MetaPhlAn2序列+索引 | 1.5G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/humann2/metaphlan2.tar.gz) | 物种注释、序列和bowtie2索引(完整版列)
MetaPhlAn2序列 | 253M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/humann2/metaphlan2.tar.gz) | 物种注释、序列，第一次运行时会自动建索引，要求安装和使用用户为同一人，有权限生成bowtie2索引
HUMAnN2 Chocophlan v0.1.1 | 5.4G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/humann2/chocophlan/full_chocophlan_plus_viral.v0.1.1.tar.gz) | 微生物泛基因组(nucleotide)，建立功能与物种组成的联系
HUMAnN2 Mapping 1.1 | 621M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/humann2/utility_mapping/full_mapping_1_1.tar.gz) | 功能描述(utility_mapping)
HUMAnN2 Uniref90 1.1 | 6.2G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/humann2/uniref/uniref90_annotated_1_1.tar.gz) | UniRef蛋白(protein)序列90%相似度聚类diamond索引
HUMAnN2 Uniref50 1.1 | 2.7G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/humann2/uniref/uniref90_annotated_1_1.tar.gz) | UniRef蛋白(protein)序列50%相似度聚类diamond索引

- MetaPhlAn4软件+数据库

主页：https://github.com/biobakery/MetaPhlAn

简介：宏基因组数据有参快速物种定量软件

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
MetaPhlAn4 | 371M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/metaphlan4.tar.gz) | metaphlan4.0.6
MetaPhlAn4 vJan21数据库 | 13G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/metaphlan4/mpa_vJan21_CHOCOPhlAnSGB_202103_bt2.tar.gz) | 物种注释、序列和bowtie2索引(完整版列)
MetaPhlAn4 vOct22数据库 | 14G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/metaphlan4/mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar.gz) | 物种注释、序列和bowtie2索引(完整版列)

### LEfSe生物标志物鉴定

主页：https://github.com/SegataLab/lefse

简介：生物标志物鉴定，结果可绘制为柱状图、物种树等。https://bioconda.github.io/recipes/lefse/README.html

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
LEfSe 1.1.2 | 441M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/lefse.tar.gz) | 生物标志物鉴定，结果可绘制为柱状图、物种树等

### Kraken2物种注释

主页：https://github.com/DerrickWood/kraken2

简介：快速的宏基因组物种分类注释软件。内存需求较大，最小16GB，推荐256GB。标准库包括人类、细菌、古菌、病毒和载体索引。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
Kraken2.1.2 | 527M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/kraken2.tar.gz) | 打包的kraken2流程，包括kraken2、bracken、krakentools、krona、r-optparse等
迷你库PlusPF | 16G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/kraken2/k2_pluspf_16gb_20230314.tar.gz) | 标准库+原生生物+真菌，压缩为16G，注释效率仅为85%，适合内存不到百G的服务器
PlusPF | 69G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/kraken2/k2_pluspf_20230314.tar.gz) | 标准库+原生生物+真菌
PlusPFP | 144G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/kraken2/k2_pluspfp_20230314.tar.gz) | 标准库+原生生物+真菌+植物

### 宏基因组组装megahit/spades

简介：宏基因组组装流程，包括megahit v1.2.9、metaSPAdes v3.15.4、MetaQUAST v5.0.2、CD-HIT v4.8.1、EMBOSS v6.6、salmon v1.8等软件。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
megahit/spades组装流程 | 768M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/megahit.tar.gz) | 宏基因组组装流程，包括megahit v1.2.9、metaSPAdes v3.15.4、MetaQUAST v5.0.2、CD-HIT v4.8.1、EMBOSS v6.6、salmon v1.8等软件



### 功能注释eggNOG/CAZy

- eggNOG综合功能注释KO/COG/CAZy/GO等

主页：http://eggnog-mapper.embl.de/

简介：免费全面的功能注释软件和数据库，可实现KO/COG/CAZy/GO等10余种注释信息，使用方便。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
eggnog-mapper 2.1.10 | 227M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/eggnog.tar.gz) | 免费全面的功能注释软件，可实现KO/COG/CAZy/GO等10余种注释信息，需要配合数据库使用，KEGG注释的替换工具
eggnog蛋白索引 | 12G | 2023-04-04 | [dmnd.gz](ftp://download.nmdc.cn/tools/meta/eggnog/eggnog.tar.gz) | eggnog蛋白注释数据库

- CAZy碳水化合物活性酶数据库

主页：http://www.cazy.org/ 

简介：最全面的碳水化合物活性酶注释数据库，扩展分析网站有dbCAN2 http://bcb.unl.edu/dbCAN2

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
dbCAN2-CAZyDB | 693M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/dbcan2/CAZyDB.08062022.tar.gz) | dbCAN2发布的整理CAZy数据库，每年更新1次


### 耐药基因 CARD

主页：https://card.mcmaster.ca/

RGI Github: https://github.com/arpcard/rgi

简介：耐药基因数据库 CARD（Comprehensive Antibiotic Resistance Database），提供了与抗菌素耐药性的分子基础相关的数据、模型和算法。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
RGI 5.2.1 | 329M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/rgi.tar.gz) | 抗生素抗性/耐药基因注释软件
Database | 3.9M| 2023-04-04 | [tar.bz2](ftp://download.nmdc.cn/tools/meta/card/data) | 抗生素抗性/耐药基因数据库


### MetaWRAP分箱

主页：https://github.com/bxlab/metaWRAP

简介：最好用的分箱流程，包括组装、提纯、定量、物种和功能注释全流程。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
MetaWRAP 1.3 | 1.7G | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/metawrap.tar.gz) | 分箱流程，包括组装、提纯、定量
MetaWRAP测试数据 | 3G | 2023-04-04 | [*.gz](ftp://download.nmdc.cn/tools/meta/metawrap) | 组装结果final.contigs.fa.gz和三个样本ERR011347/8/9的双端文件_1/2.fastq.gz
CheckM | 288M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/checkm/checkm_data_2015_01_16.tar.gz) | 基因组完整性评估，用于binning, bin_refinement, reassemble_bins
NCBI_nt | 201GB | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/NCBI/nt) | 基因组物种注释，用于blobology, classify_bins
NCBI_tax | 58M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/NCBI/tax/taxdump.tar.gz) | 基因组物种注释，用于blobology, classify_bins

### dRep基因组去冗余

GitHub: https://github.com/MrOlm/drep

Conda: https://bioconda.github.io/recipes/drep/README.html

简介：细菌基因组、宏基因组组装基因组去冗余，鉴定种、株级别的基因组数量。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
dRep v3.2.2 | 246M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/drep.tar.gz) | 细菌基因组、宏基因组组装基因组去冗余
CheckM | 288M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/checkm/checkm_data_2015_01_16.tar.gz) | 基因组完整性评估，用于drep，同metawrap中CheckM


### GTDB细菌基因组物种注释

网址：https://gtdb.ecogenomic.org/

该数据库于2018/2020连发两篇Nature Biotechnology，被引数千次，是引领全基因组细菌、古菌分类中最新的数据库、最前沿的方法。

版本 | 大小 | 更新时间 | 下载链接 | 说明 
-|-|-|-|-
GTDB-tk v2.2.6 | 154M | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/conda/gtdbtk.tar.gz) | 细菌、古菌基因组物种注释
GTDB r207v2 | 67GB | 2023-04-04 | [tar.gz](ftp://download.nmdc.cn/tools/meta/gtdb/gtdbtk_r207_v2_data.tar.gz) | 分箱、细菌基因组物种注释、进化树构建