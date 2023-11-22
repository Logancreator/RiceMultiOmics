library(methylKit)
library(dplyr)
library(optparse)
 
#step1 设置参数格式及帮助文档 make_option
# 描述参数的解析方式
option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = FALSE,
              action = "store", help = "This is file!"
  ),
  make_option(c("-n", "--name"), type = "character", default = FALSE,
              action = "store", help = "This is name!"
  )
)

 
#type：指定参数的数据类型，可以是 “logical”, “integer”, “double”, “complex”, or “character”
#default：指定选项的默认值
#action：指定拿到参数后的动作，可以是store（存储为字符），store_true（存储为逻辑值True），store_false（存储为逻辑值False）
#help：指定帮助信息
 
#step2 构建参数对象 OptionParser,以便解析
opt_parser = OptionParser(
      usage = "usage: %prog [options]",
      option_list = option_list,
      add_help_option = TRUE,
      prog=NULL , 
      description = "This Script is to calculate statistics value.")
 
#usage：打印使用说明
#option_list：输入第一步构建的参数列表
#description:打印的帮助信息。
#add_help_option :是否自动添加 -h/--help 参数，默认添加
#prog：指明脚本名称，可以不填写，当在命令行中使用时，程序会自动识别脚本名称。
#description:在options语句之后打印的print_help的附加文本 
 
 
#step3 解析参数 parse_args
opt = parse_args(opt_parser);
 
#step4 提取参数的值

file.list=list(opt$file)
name=list(opt$name)


# read the files to a methylRawList object: myobj
myobj=processBismarkAln(location=file.list,
                       sample.id=name,assembly="goose",
                       save.folder="/public/home/changjianye/project/shenxu_data/geese_MBH_WGBS/Batmeth2/BatMeth3_result/tmp",
                       save.context=c("CpG","CHG","CHH"),read.context=c("CpG","CHG","CHH"),
                       save.db=TRUE,nolap=FALSE,mincov=5,minqual=20,treatment=c(1))