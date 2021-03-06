gtts\_de\_limma（limma差异表达分析工具）
----------------------------------------

### 介绍

limma是一个常用的差异表达分析工具，主要用于芯片表达谱数据。本程序通过limma对芯片表达谱进行差异表达分析。

### 开发环境

1.  系统版本：macOS catalina 10.15.3
2.  R版本：3.6.2

### 依赖

1.  argparser: 0.6
2.  limma: 3.42.2
3.  logger: 0.1

### 输出信息

#### 警告信息

程序运行过程中，由于用户的输入错误导致的警告，已定义的有：

    kWarns <- list(
      WARN_GTTS_DE_LIMMA_001 = '没有差异基因'
    )

#### 错误说明

程序运行过程中，由于用户的输入错误导致的错误，已定义的有：

    kErrs <- list(
      ERR_GTTS_DE_LIMMA_001 = '未指定表达谱文件路径',
      ERR_GTTS_DE_LIMMA_002 = '表达谱文件不存在',
      ERR_GTTS_DE_LIMMA_003 = '表达谱文件格式不支持',
      ERR_GTTS_DE_LIMMA_004 = '表达谱文件格式未能正确读取',
      ERR_GTTS_DE_LIMMA_005 = '未指定表型文件路径',
      ERR_GTTS_DE_LIMMA_006 = '表型文件不存在',
      ERR_GTTS_DE_LIMMA_007 = '表型信息文件格式不支持',
      ERR_GTTS_DE_LIMMA_008 = '表型信息文件格式未能正确读取',
      ERR_GTTS_DE_LIMMA_009 = 'log fold change必须大于0',
      ERR_GTTS_DE_LIMMA_010 = 'p value必须大于等于0小于等于1',
      ERR_GTTS_DE_LIMMA_011 = '非法校正方法',
      ERR_GTTS_DE_LIMMA_012 = '输出文件夹不存在',
      ERR_GTTS_DE_LIMMA_013 = '输出文件夹权限错误'
    )

### 参数

    usage: gtts_de_limma.R [--] [--help] [--debug] [--opts OPTS]
           [--exp-file EXP-FILE] [--exp-file-fmt EXP-FILE-FMT]
           [--pheno-file PHENO-FILE] [--pheno-file-fmt PHENO-FILE-FMT]
           [--grp-col-num GRP-COL-NUM] [--lfc LFC] [--pvalue PVALUE]
           [--adjust-method ADJUST-METHOD] [--outdir OUTDIR]

    limma差异表达分析

    flags:
      -h, --help           show this help message and exit
      -d, --debug          打开调试模式

    optional arguments:
      -x, --opts           RDS file containing argument values
      -e, --exp-file       表达谱文件路径
      --exp-file-fmt       表达谱文件格式, 支持tsv,csv,rds [default: tsv]
      -p, --pheno-file     临床信息文件路径
      --pheno-file-fmt     表型信息文件格式, 支持tsv,csv,rds [default: tsv]
      -g, --grp-col-num    表型信息文件中哪一列是样本分组信息 [default: 2]
      -l, --lfc            log fold change [default: 1]
      --pvalue             p value [default: 0.05]
      -a, --adjust-method  校正方法，支持none,BH,BY,holm [default: BH]
      -o, --outdir         输出文件夹 [default: .]

1.  `--exp-file`

    表达谱文件路径：

         1. 不能为空
         2. 行是基因，且必须有行名
         3. 列是样本，且必须有列名

2.  `--exp-file-fmt`

    表达谱文件格式：

         1. 可以为空，当为空时默认值为tsv格式
         2. 可选的格式有tsv、csv、rds

3.  `--pheno-file`

    表型信息文件路径：

         1. 不能为空
         2. 行是样本，且必须有行名
         3. 列是样本表型信息项，且必须有列名

4.  `--pheno-file-fmt`

    表型信息文件格式：

         1. 可以为空，当为空时默认值为tsv格式
         2. 可选的格式有tsv、csv、rds

5.  `--grp-col-num`

    表型信息文件中哪一列是分组信息（行名计算在内）：

         1. 可以为空，当为空时默认值为2
         2. 因为第一列为行名，所以至少为第2列

6.  `--lfc`

    log fold change值的阈值，大于该值的基因会被输出：

         1. 可以为空，当为空时，默认值为1
         1. 必须大于0

7.  `--pvalue`

    p值的阈值，小于该值的基因会被输出：

         1. 可以为空，当为空时，默认值为0.05

8.  `--adjust-method`

    校正的方法：

         1. 可以为空，当为空时，默认值为BH
         2. 支持的方法有：
             1. none，即不校正
             2. BH
             3. BY
             4. holm

9.  `--outdir`

    结果的输出路径，当有多个分组时，会有多个结果输出，所以必须指定路径：

         1. 可以为空，当为空时，默认值为`.`，即当前文件夹

### 结果说明

假如，表型信息中的分组有三个分组（case\_1、case\_2、control），那么会产生三个文件：

    1. case_1-case_2.tsv：case_1和case_2之间的差异基因
    2. case_1-control.tsv：case_1和control之间的差异基因
    3. case_2-control.tsv：case_2和control之间的差异基因

以`case_1-case_2.tsv`为例，该文件4列：

    ""      "logFC" "P.Value"       "adj.P.Val"
    "TIGD1" 1.05399384016666        0.00315109772350737     0.977683104849865
    "ENTPD3-AS1"    -0.789166808000002      0.00394961089028577     0.977683104849865
    "LOC284578"     1.071551847     0.00405965900488913     0.977683104849865
    "LINC00847"     0.784985876333331       0.005284773539729       0.977683104849865
    "LINC00317"     -0.285425680166669      0.00680254636644037     0.992011277422333
    "CTA-390C10.10" -0.741394700000003      0.0102304714915927      0.992011277422333
    "MBLAC1"        0.541562605833329       0.0149161401177375      0.992011277422333
    "PDILT" 0.378897207499999       0.015206389960181       0.992011277422333
    "TAAR9" 0.350609278333332       0.0152704495850186      0.992011277422333

1.  第一列为差异基因
2.  第二列为log fold change值
3.  第三列为校正前的p值
4.  第四列时校正后的p值
