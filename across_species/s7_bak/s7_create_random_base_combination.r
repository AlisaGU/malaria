#!/home/gushanshan/anaconda3/envs/vscode_r/bin/Rscript
# 1. packages and external scripts ---------------------------------------- TODO:


# 2. functions ------------------------------------------------------------ TODO:


# 3. input ---------------------------------------------------------------- TODO:
Args <- commandArgs()
whole_region <- as.integer(Args[6])
A_count <- as.integer(Args[7])
T_count <- as.integer(Args[8])
C_count <- as.integer(Args[9])
G_count <- as.integer(Args[10])

# 4. variable setting of test module--------------------------------------- TODO:


# 5. process -------------------------------------------------------------- TODO:
A_pro <- A_count / whole_region
T_pro <- T_count / whole_region
C_pro <- C_count / whole_region
G_pro <- G_count / whole_region
set.seed(1)
seq <- sample(c("A", "T", "C", "G"), 6, replace = T, prob = c(A_pro, T_pro, C_pro, G_pro))
print(seq)
