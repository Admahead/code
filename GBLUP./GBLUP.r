###GBLUP贝叶斯模型

# 加载必要的包
library(BGLR)
library(genio)
install.packages("genio")

# 1. 读取PLINK文件
plink_data <- read_plink("geno")  # 不需要后缀，会自动读取.bed/.bim/.fam

# 提取基因型矩阵 (个体 × SNP，编码为0,1,2)
geno_matrix1 <- plink_data$X
dim(geno_matrix1)
geno_matrix2 <- t(geno_matrix1)
 
# 2. 读取表型数据
pheno <- read.table("phe.txt", header = TRUE)
names(pheno) <- c("ID", "Gender", "Birth_Date", "Test_Farm", 
                  "Daily_Gain", "Daily_Feed_Intake", "Feed_Conversion_Ratio")

# 3.计算基因组关系矩阵（G矩阵）- GBLUP的核心


p2 <- colMeans(geno_matrix2)/2
scale_factor2 <- sum(2*p2*(1-p2))
M2 <- scale(geno_matrix2,center = TRUE,scale = FALSE)
G_matrix <- tcrossprod(M2) / scale_factor2

# 设置GBLUP模型参数
nIter <- 12000    # 总迭代次数
burnIn <- 200    # 预烧期

# 运行GBLUP模型
# 使用RKHS（Reproducing Kernel Hilbert Space）方法实现GBLUP
ETA <- list(list(K = G_matrix, model = "RKHS"))

# 运行模型
set.seed(123)  # 设置随机种子保证结果可重复
fm_gblup <- BGLR(y = pheno$Feed_Conversion_Ratio, 
                 ETA = ETA, 
                 nIter = nIter, 
                 burnIn = burnIn, 
                 saveAt = "gblup_",
                 verbose = TRUE)

# 提取结果和计算遗传参数# 提取方差组分
Vg <- fm_gblup$ETA[[1]]$varU  # 遗传方差
Ve <- fm_gblup$varE           # 残差方差
Vp <- Vg + Ve                 # 表型方差

# 计算基因组遗传力
h2 <- Vg / Vp



#读文件
geno_raw <- read.table("GENO.raw", header = TRUE, sep = "")
geno_matrix <- as.matrix(geno_raw[, -c(1:6)])
dim(geno_matrix)



p <- colMeans(geno_matrix)/2
scale_factor <- sum(2*p*(1-p))
M <- scale(geno_matrix,center = TRUE,scale = FALSE)


# 设置BayesA模型参数
nIter <- 12000    # 总迭代次数
burnIn <- 2000    # 预烧期


# 运行模型BayesA
set.seed(123)  # 设置随机种子保证结果可重复
ETA_A <- list(list(X = M-1 , model ="BayesA" ))
fm_bayesA <- BGLR(y = pheno$Feed_Conversion_Ratio, 
                 ETA = ETA_A, 
                 nIter = nIter, 
                 burnIn = burnIn, 
                 saveAt = "bayesA_",
                 verbose = TRUE)

Ve_A <- fm_bayesA$varE
Vg_A <- sum(fm_bayesA$ETA[[1]]$b ^2)
Vp_A <- Ve_A +Vg_A
h2_A <- Vg_A /Vp_A


