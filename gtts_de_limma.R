#!/usr/bin/env Rscript
# copyright(c) 2020 gaptech
# author：gttlixing
# 使用limma对芯片数据进行差异表达分析
library(argparser)
library(limma)
library(logger)


# 警告类型预定义
kWarns <- list(
  WARN_GTTS_DE_LIMMA_001 = '没有差异基因'
)
# 错误类型预定义
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
# 表达谱支持的文件格式
kSupportedExpFileFmts <- c(
  'tsv', 'csv', 'rds'
)
# 默认表达谱文件格式
kDefaultExpFileFmt <- 'tsv'
# 表型信息文件支持的文件格式
kSupportedPhenoFileFmts <- c(
  'tsv', 'csv', 'rds'
)
# 默认表达谱文件格式
kDefaultPhenoFileFmt <- 'tsv'
# 分组字符如果含有空格的替换值
kReplaceSpaceWith <- '_'
# 多个分组之间的间隔符
kDelimiterOfGrps <- '-'
# 默认哪一列是分组信息
kDefaultGrpColNum <- 2
# 默认log fold change
kDefaultLfc <- 1
# 默认p value
kDefaultPvalue <- 0.05
# 校正方法
kAdjustMethods <- c(
  'none', 'BH', 'BY', 'holm'
)
kDefaultAdjustMethod <- 'BH'
# 结果中保留哪些列
kKeepCols <- c(
  'logFC', 'P.Value', 'adj.P.Val'
)
# 默认输出路径
kDefaultOutdir <- '.'


GetMsg <- function(code) {
  # 获取格式化的警告或错误信息
  # 一般一个格式化的警告错误信息由两个部分组成
  # 警告或错误代号：警告或错误信息
  # Args:
  #   code: 警告或错误代号
  # Returns:
  #   格式化的警告或错误信息
  type <- substr(code, 1, 3)
  if (type == 'ERR') {
    msgs <- kErrs
  } else {
    msgs <- kWarns
  }
  return(sprintf('%s：%s', code, msgs[[code]]))
}


ReadExpFile <- function(file, fmt = 'tsv') {
  # 读取表达谱
  # Args:
  #   file: 表达谱文件路径
  #   fmt: 表达谱文件格式，支持格式tsv,csv,rds
  # Returns:
  #   表达谱数据框
  exp <- NULL
  if (fmt == 'tsv') {
    exp <- read.delim(
      file, row.names = 1,
      check.names = FALSE
    )
  } else if (fmt == 'csv') {
    exp <- read.csv(
      file, row.names = 1,
      check.names = FALSE
    )
  } else if (fmt == 'rds') {
    exp <- readRDS(file)
  }
  if (is.null(exp)) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_004'))
  }
  return(exp)
}


ReadPhenoFile <- function(file, fmt = 'tsv') {
  # 读取表型信息文件
  # Args:
  #   file: 表型信息文件路径
  #   fmt: 表型信息文件格式，支持tsv,csv,rds
  # Returns:
  #   样本表型信息数据框
  pheno <- NULL
  if (fmt == 'tsv') {
    pheno <- read.delim(
      file, row.names = 1,
      check.names = FALSE
    )
  } else if (fmt == 'csv') {
    pheno <- read.csv(
      file, row.names = 1,
      check.names = FALSE
    )
  } else if (fmt == 'rds') {
    pheno <- readRDS(file)
  }
  if (is.null(pheno)) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_008'))
  }
  return(pheno)
}


MakeDesign <- function(grp) {
  # 获取design
  # Args:
  #   grp: 命名分组因子
  # Returns:
  #   design矩阵
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  rownames(design) <- names(grp)
  return(design)
}


MakeContrastsMatrix <- function(grp) {
  # 创建contrast matrix 
  # Args:
  #   grp: 命名分组因子
  # Returns:
  #   contrast matrix
  grp.combn <- combn(levels(grp), 2)
  contrasts <- paste0(grp.combn[1,], kDelimiterOfGrps, grp.combn[2,])
  return(makeContrasts(
    contrasts = contrasts,
    levels = levels(grp)
  ))
}


DoDe <- function(exp, design, contrastsMatrix) {
  # 差异表达分析
  # Args:
  #   exp: 表达谱数据框
  #   design: design
  #   contrastsMatrix: contrasts matrix
  # Returns:
  #   bayes fit result
  vfit <- lmFit(exp, design)
  vfit <- contrasts.fit(vfit, contrasts = contrastsMatrix)
  return(eBayes(vfit))
}


GetDegs <- function(efit, lfc = 1, pValue = 0.05, adjustMethod = 'BH') {
  # 获取差异表达基因
  # Args:
  # Returns:
  dt <- decideTests(efit, lfc = lfc, p.value = pValue,
                    adjust.method = adjustMethod)
  apply(dt, 1, function(gene) any(gene != 0))
  degs <- names(which(apply(dt, 1, function(gene) any(gene != 0))))
}


GetArgs <- function() {
  # 获取命名行参数
  # Args:
  # Returns:
  #   命令行参数list
  parser <- arg_parser('limma差异表达分析')
  parser <- add_argument(parser, '--exp-file',
                         help = '表达谱文件路径')
  parser <- add_argument(parser, '--exp-file-fmt',
                         default = kDefaultExpFileFmt,
                         help = sprintf(
                           '表达谱文件格式, 支持%s',
                           paste(kSupportedExpFileFmts, collapse = ',')
                         ))
  parser <- add_argument(parser, '--pheno-file',
                         help = '临床信息文件路径')
  parser <- add_argument(parser, '--pheno-file-fmt',
                         default = kDefaultPhenoFileFmt,
                         help = sprintf(
                           '表型信息文件格式, 支持%s',
                           paste(kSupportedPhenoFileFmts, collapse = ',')
                         ))
  parser <- add_argument(parser, '--grp-col-num',
                         default = kDefaultGrpColNum,
                         help = '表型信息文件中哪一列是样本分组信息')
  parser <- add_argument(parser, '--lfc',
                         default = kDefaultLfc,
                         help = 'log fold change')
  parser <- add_argument(parser, '--pvalue',
                         default = kDefaultPvalue,
                         help = 'p value')
  parser <- add_argument(parser, '--adjust-method',
                         default = kDefaultAdjustMethod,
                         help = sprintf(
                           '校正方法，支持%s',
                           paste(kAdjustMethods, collapse = ',')
                         ))
  parser <- add_argument(parser, '--outdir', default = kDefaultOutdir,
                         help = '输出文件夹')
  parser <- add_argument(parser, '--debug', flag = TRUE, help = '打开调试模式')
  return(parse_args(parser))
}


CheckArgs <- function(args) {
  # 命令行参数校验
  # 检查是否指定表达谱文件路径
  if (is.na(args$exp_file)) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_001'))
  }
  
  # 检查表达谱文件是否存在
  if (! file.exists(args$exp_file)) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_002'))
  }
  
  # 检查表达谱文件格式是否支持
  if (! args$exp_file_fmt %in% kSupportedExpFileFmts) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_003'))
  }
  
  # 检查是否指定表型信息文件路径
  if (is.na(args$pheno_file)) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_005'))
  }
  
  # 检查表型文件是否存在
  if (! file.exists(args$pheno_file)) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_006'))
  }
  
  # 检查表型文件格式是否支持
  if (! args$pheno_file_fmt %in% kSupportedPhenoFileFmts) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_007'))
  }
  
  # 检查lfc是否大于0
  if (args$lfc < 0) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_009'))
  }
  
  # 检查p value是否大于等于0小于等于1
  if (args$pvalue < 0 | args$pvalue > 1) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_010'))
  }
  
  # 检查adjust method是否支持
  if (! args$adjust_method %in% kAdjustMethods) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_011'))
  }
  
  outdir <- normalizePath(args$outdir)
  # 检查输出文件夹是否存在
  if (! dir.exists(outdir)) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_012'))
  }
  
  # 检查输出文件夹是否可写入
  if (file.access(outdir, 2) != 0) {
    stop(GetMsg('ERR_GTTS_DE_LIMMA_013'))
  }
  return(args)
}


Main <- function() {
  # 脚本入口函数
  args <- GetArgs()
  args <- CheckArgs(args)
  if (args$debug) {
    log_threshold(DEBUG)
  } else {
    log_threshold(WARN)
  }
  outdir <- normalizePath(args$outdir)
  grp.col.num <-  args$grp_col_num - 1  # 第一列是列名

  # 读取表达谱和表型信息
  exp <- ReadExpFile(args$exp_file, args$exp_file_fmt)
  pheno <- ReadPhenoFile(args$pheno_file, args$pheno_file_fmt)

  # 取交集并获取分组
  log_debug(sprintf('表达谱文件中共有样本%s个', dim(exp)[2]))
  log_debug(sprintf('表型信息文件中共有样本%s个', dim(pheno)[1]))
  samples <- intersect(colnames(exp), rownames(pheno))
  log_debug(sprintf('既有表达谱信息、又有表型信息的样本共有%s个', length(samples)))
  exp <- exp[, samples]
  grp <- as.character(pheno[samples, grp.col.num])
  grp <- make.names(grp)  # 防止分组里面有非法字符
  names(grp) <- samples
  grp <- as.factor(grp)

  # 差异表达分析
  design <- MakeDesign(grp)
  contrasts.matrix <- MakeContrastsMatrix(grp)
  efit <- DoDe(exp, design, contrasts.matrix)

  # 获取差异表达基因并保存到文件
  degs <- lapply(colnames(contrasts.matrix), function(coef) {
    tt <- topTable(
      efit,
      coef = coef,
      number = Inf,
      adjust.method = args$adjust_method,
      p.value = args$pvalue,
      lfc = args$lfc
    )
    if (dim(tt)[1] < 1) {  # 没有差异基因显示警告
      log_warn(GetMsg('WARN_GTTS_DE_LIMMA_001'))
    } else {
      write.table(
        tt[kKeepCols],
        file.path(outdir, sprintf('%s.tsv', coef)),
        col.names = NA,
        sep = '\t'
      )
    }
    return(tt)
  })
}

Main()
