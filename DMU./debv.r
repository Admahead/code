# =========================
# 1. 加载所需R包
# =========================
# blupADC：用于家系处理、调用DMU等遗传评估相关分析
# dplyr：用于数据整理、筛选、合并等操作
library(blupADC)
library(dplyr)


# =========================
# 2. 设置工作目录和结果输出目录
# =========================
# 获取当前工作目录
base_dir <- getwd()

# 在当前工作目录下新建一个结果文件夹，用来存放DMU输出结果
output_result_path <- file.path(base_dir, "PBLUP_GLLBody")

# 如果文件夹不存在则创建；recursive = TRUE表示可递归创建多级目录
dir.create(output_result_path, recursive = TRUE, showWarnings = FALSE)


# =========================
# 3. 读取系谱文件（Pedigree）
# =========================
# 读取LL群体的系谱文件
# header = TRUE 表示第一行为列名
# stringsAsFactors = FALSE 防止字符型变量自动转成因子
Ped <- read.table(
  "/data6/zlai/4-OtherTrait.zip/4-OtherTrait/4-OtherTrait/2-Trait/DEBV_LL/LL.ped",
  header = TRUE,
  stringsAsFactors = FALSE
)


# =========================
# 4. 整理系谱并生成DMU可用格式
# =========================
# trace_pedigree() 用于追踪和整理系谱关系
# output_pedigree_tree = TRUE 表示输出完整家系树结构
pedigree_result <- trace_pedigree(
  input_pedigree = Ped,
  output_pedigree_tree = TRUE
)

# 提取重编码后的系谱表
# 一般会把原始ID转换为DMU更适合使用的数字ID
Ped_rename <- pedigree_result$rename_ped

# 保存重编码后的系谱表，便于后续核对
write.table(
  Ped_rename, "Ped.rename.txt",
  row.names = FALSE, quote = FALSE, sep = "\t"
)


# =========================
# 5. 构建DMU所需的系谱关系文件
# =========================
# 从重编码后的系谱表中提取第3到第6列
# 这里通常对应：个体ID、父本ID、母本ID、出生顺序/群体等信息
Ped_Ana <- Ped_rename[, 3:6]

# 新增一列 Order，赋值为0
# 这个通常是DMU分析需要的附加列
Ped_Ana$Order <- 0

# 保存为DMU分析使用的系谱文件
# col.names = FALSE 表示不写列名，因为DMU通常要求纯数据格式
write.table(
  Ped_Ana, "Ped_Ana.txt",
  row.names = FALSE, col.names = FALSE,
  quote = FALSE, sep = "\t"
)


# =========================
# 6. 建立原始ID与数字ID的对应表
# =========================
# 从重编码后的系谱文件中提取两列：
# 第1列：原始个体ID
# 第3列：转换后的数字ID
ped_map <- Ped_rename[, c(1, 3)]

# 修改列名，方便后面和表型文件做匹配
colnames(ped_map) <- c("original_id", "numeric_id")


# =========================
# 7. 读取表型文件
# =========================
# 读取LL群体的BodyL性状表型文件
phe0 <- read.table(
  "BodyL.txt",#表型文件，包括固定效应等
  header = TRUE,
  stringsAsFactors = FALSE
)


# =========================
# 8. 表型数据与系谱数字ID进行匹配
# =========================
# 用左连接方式把表型文件中的ID和ped_map中的original_id对应起来
# 得到每个个体对应的numeric_id
phe1 <- phe0 %>%
  left_join(ped_map, by = c("ID" = "original_id")) %>%
  
  # 去掉没有匹配上数字ID的个体
  filter(!is.na(numeric_id)) %>%
  
  # 去掉BodyL缺失的个体
  filter(!is.na(BodyL))


# =========================
# 9. 构建DMU分析所需的表型文件
# =========================
# 将原始ID替换为numeric_id，并只保留分析需要的列
phe_dmu <- phe1 %>%
  mutate(ID = numeric_id) %>%
  select(ID, BFactory, Byear, Bseason, Sex, Parity, BodyL)

# 保存DMU使用的表型文件
# col.names = FALSE 表示不输出列名，符合DMU输入格式要求
write.table(
  phe_dmu, "P_BodyL.txt",
  row.names = FALSE, col.names = FALSE,
  quote = FALSE, sep = "\t"
)

# 另外保存一份列名文件，方便自己核对每一列代表什么
write.table(
  colnames(phe_dmu), "names_BodyL.txt",
  row.names = FALSE, col.names = FALSE,
  quote = FALSE
)


# =========================
# 10. 调用DMU进行PBLUP_A分析
# =========================
# run_DMU() 是blupADC中用于自动生成DMU参数文件并运行分析的函数
run_DMU(
  # 表型文件中的列名顺序
  phe_col_names = c("ID","BFactory","Byear","Bseason","Sex","Parity","BodyL"),

  # 目标性状名称
  target_trait_name = list("BodyL"),

  # 固定效应：场次、年份、季节、性别、胎次
  fixed_effect_name = list(
    c("BFactory","Byear","Bseason","Sex","Parity")
  ),

  # 协变量效应：这里没有设置，所以为NULL
  covariate_effect_name = NULL,

  # 随机效应名称：这里设置个体ID为随机遗传效应
  random_effect_name = list(c("ID")),

  # 指定遗传效应对应的列名
  genetic_effect_name = "ID",

  # 是否加入永久环境效应
  # 这里为FALSE，表示不加入永久环境效应
  included_permanent_effect = list(FALSE),

  # 表型文件名称
  phe_name = "P_BodyL.txt",

  # 表型文件所在路径
  phe_path = base_dir,

  # integer_n = 6 表示前6列按整数型处理
  # 通常是ID和分类固定效应列
  integer_n = 6,

  # 分析模型：PBLUP_A 表示基于系谱关系矩阵A的PBLUP
  analysis_model = "PBLUP_A",

  # DMU模块：dmuai 常用于方差组分和育种值估计
  dmu_module = "dmuai",

  # 系谱关系文件名称
  relationship_name = "Ped_Ana.txt",

  # 系谱关系文件所在路径
  relationship_path = base_dir,

  # 是否计算dEBV
  cal_debv = TRUE,

  # 计算dEBV时所用的系谱文件
  debv_pedigree_name = "Ped_Ana.txt",
  debv_pedigree_path = base_dir,

  # 指定输出结果目录
  output_result_path = output_result_path
)
