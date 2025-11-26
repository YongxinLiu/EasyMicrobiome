#!/usr/bin/perl -w
# 加载时间管理，参数管理，文件名和路径处理的基础包，无须安装
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Scripts usage and about.
# 程序的帮助文档，良好的描述是程序重用和共享的基础，也是程序升级和更新的前提
###############################################################################
sub usage {
    die(
        qq!
Usage:    format_drep2cluster.pl -i temp/drep95/data_tables/Cdb.csv -d temp/drep95/data_tables/id -o temp/drep95/data_tables/Cdb.list -h header num
Function: Format drep result into rep and list
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, liuyongxin_bio\@163.com, QQ:42789409
Version:  v1.01
Update:   2024/4/18
Notes:    v1.00 2024/4/18 Inital the script
Notes:    v1.01 2025/11/18 Update the script
\n!
    )
}

###############################################################################
#命令行参数据的定义和获取，记录程序初始时间，设置参数默认值
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:e:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
# 调置参数的初始值，可以添加更多参数的默认值
$opts{h}=0 unless defined($opts{h});

###############################################################################
#读入的数据或注释文件，用于与输入文件比较或注释(可选)，提供三种方式
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
# 1. 散列结构数据库，要求数据文件有唯一ID并且无顺序要求
my %database; #database in hash
#X026b2
#X033
#X046b2
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[1]}=1;
}
# 2. 数组结构数据库，无唯一ID，但有顺序要求
#my (@tmp1,@tmp2); #database in array
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;
# 3. 批量数据文件，读取一批有相似结构的文件
#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
# 正文部分，读取输入文件，列出输入和输入文件的三行作为示例，方便编程处理数据
open INPUT,"<$opts{i}";
#Y828.fa,1_1,0.050000000000000044,average,ANImf,1
#Y954.fa,1_1,0.050000000000000044,average,ANImf,1
#Y011b1.fa,2_1,0.050000000000000044,average,ANImf,2
#Y091b1.fa,2_1,0.050000000000000044,average,ANImf,2
#Y136b1.fa,3_0,0.050000000000000044,average,ANImf,3
#Y593.fa,4_0,0.050000000000000044,average,ANImf,4
open OUTPUT,">$opts{o}";
#X026b2	Y401b6,Y454b6,X021b4,X026b2,X029b3,X040b2,X043b3,X051b3,X092,X292b4,X293b1,X293b4,X297b2,X298b3,X299b1,X299b6,X302b1,X302b3,X330b1,X330b3,Y003b5,Y085,Y1015,Y1018,Y134,Y136b2,Y745,Y746,Y769,Y877,Y955
#X033	X015,X016b1,X018b2,X020b1,X021b1,X022b3,X024b2,X025b1,X026b1,X028b1,X029b1,X030b2,X031b1,X032b1,X033,X036b1,X037,X038,X039b1,X040b1,X041b1,X043b1,X046b1,X047b1,X048b1,X049b2,X051b1,X054b2,X063b1,X064b1,X070b4,X073b2,X074b2,X076b2,X077b2,X078b3,X080b2,X081b1,X082b1,X086b1,X359b1,X360,Y584b1,Y584b7,Y753,Y859,Y862,Y943,Y948,Y957b2,Y976,Y980b2,Y982,Y983b1
#X046b2	X224b1,X046b2,X048b3,X049b4,X059b1,X063b6,X064b4,X080b4,X303b5
#my %count;
## h参数用于去除有文件头的行
#while ($opts{h}>0) { #filter header
#	$tmp=<INPUT>;
#	$opts{h}--;
#	# 可选，输出文件也保留文件头
#	#print OUTPUT $tmp;
#}
# 输入和输入处理部分，常用按行读取处理并输入，默认按tab分割数据

# 使用存储结果，防止输出重复结果
print OUTPUT "ID\tList\n";
while (<INPUT>) {
	chomp;
	my @tmp=split/,/;
	my @tmp1=split('\.',$tmp[0]);
	print $tmp1[0],"\n";
# 新版结果没有|，跳过此步
# if ($tmp[1]=~/\|$/) {
   $i=$#tmp1;
	print $i,"\n";
# }else{
#   $i=$#tmp1-1;
# }
	foreach (1..$i) {
		# 跳过酶学编号的结果
		if ($tmp1[$_]=~/^\d/) {
			next;
		}
		# 去除下划线后面内容
		# $tmp1[$_]=~s/_.*//;
		if (!$database{$tmp[0]}{$tmp1[$_]}) {
		#print "$tmp[0]\t$tmp1[$_]\n";
		if ($tmp[10] < $opts{e}) {
  		print OUTPUT "$tmp[0]\t$tmp1[$_]\t$tmp[2]\t$tmp[10]\n";
		}
		}else{
			next;
		}
		$database{$tmp[0]}{$tmp1[$_]}++;
	}
#	print OUTPUT "$tmp[0]\t$tmp[1]\n";
}
close INPUT;


close OUTPUT;

###############################################################################
#Record the program running time!
# 输出程序运行时间
###############################################################################
my $duration_time=time-$start_time;
#print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
#print "This compute totally consumed $duration_time s\.\n";

