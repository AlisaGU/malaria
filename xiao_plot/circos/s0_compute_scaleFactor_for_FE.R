#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:


# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
### 6mA KD
## 用Pf3D7_01_v3的第327位计算
FE <- 4.41724 ## XLS表格第8列：fold_enrichment
T_cap <- 31.71773 ## treatment bedgraph文件中该位点对应的value
C_cap <- 7.08750 ## control bedgraph文件中该位点对应的value
S <- (FE - 1) / (T_cap - FE * C_cap) ## 8.323738

## 用Pf3D7_01_v3的第539位验证
T_cap <- 32.19830
C_cap <- 7.75195
FE <- (T_cap * S + 1) / (C_cap * S + 1) ## 计算结果为1.409521，与xls表格(4.10544)的一致
