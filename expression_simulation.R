# 1. Tell R which packages you need to use

library(readr)

# 2. Upload data *replace ###### with path to file "all_patient_TPM.txt"
# Before chloroquine treatment = "X_0" and After chloroquine treatment (8 hours) = "X_1"

targets <- read_delim("~/######/all_patient_TPM.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
View(targets)

# 3. Remove samples without chloroquine treatment data

RelapseDrops <- c("AR_0","AR_1","BR_0","BR_1","CR_0","CR_1","DR_0","DR_1","ER_0","ER_1","FR_0","FR_1","KR_0","KR_1","AR_1")

targets <- targets[,!(names(targets) %in% RelapseDrops)]
targets.a <- targets[,c(5:44)]

# 4. Rename rows using the gene names

row.names(targets.a) <- targets$`GeneID`

targets.Basal <- targets[,c("A_0","B_0","C_0","D_0","E_0","F_0","G_0","H_0","I_0","J_0","K_0","L_0","M_0","N_0","O_0","P_0","Q_0","S_0","T_0","U_0")]

# 5. Initialize vector for all iterations of the simulation (100 iterations total)
total = vector(,101)

# 6. Simulation: Generate random vector of samples after chloroquine treatment, compute pairwise gene comparisons, and sum

i = 1
while (i<=100) {
  y_names <- c("A_1","B_1","C_1","D_1","E_1","F_1","G_1","H_1","I_1","J_1","K_1","L_1","M_1","N_1","O_1","P_1","Q_1","S_1","T_1","U_1")
  y_rand <- sample(y_names, 20, replace=F)
  value.iter = vector(,20)
  #print(paste(value.iter))
  for (col_y in 1:ncol(targets.Basal)) {
    x <- y_rand[col_y]
    Basal <- targets.Basal[col_y]
    CQ <- targets.a[x]
    value.iter[c(col_y)] <- sum((Basal-CQ)/(Basal+1))
  }
  #print(paste(value.iter))
  total[i] <- sum(value.iter)
  test <- sum(value.iter)
  i = i+1
}

# 7. Simulation: Compute pairwise gene comparisons and sum for properly paired samples

for (col_y in 1:ncol(targets.Basal)) {
  x <- 1
  Basal <- targets.Basal[col_y]
  CQ <- targets.a[x]
  value.iter[c(col_y)] <- sum((Basal-CQ)/(Basal+1))
  x = x+1
}
total[i] <- sum(value.iter)
test <- sum(value.iter)

# 8. Write file containing values. Last entry is the correct comparison 

write.table(total, na="NA", file="expression_comparison_simulations.txt")
