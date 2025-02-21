#!/usr/bin/env Rscript

library(lme4)
library(dplyr)
library(xlsx)

# read growth data
df <- paste('../../amiga-clade-5/summary/merged_summary.txt')
df <- read.table(df,sep="\t",header=TRUE)

# to harmonize analyses, rename column headers
names(df)[names(df) == 'Carbon.Source'] <- 'Substrate'
names(df)[names(df) == 'k_lin'] <- 'Carrying_Capacity'

# define substrates to analyse
simple_sugars <- c("Glucose", "Fructose", "Tagatose", "Ribose")
other_substrates <- c("N-acetylneuraminic acid", "N-acetylglucosamine", "Mannitol", "Salicin")
substrates <- c(simple_sugars, other_substrates)

# define substrate groups
df <- df %>%
  mutate(Substrate_Group = case_when(
    Substrate %in% simple_sugars ~ "Simple_Sugar",
    Substrate %in% other_substrates ~ "Other_Substrate",
    TRUE ~ NA_character_
  ))

# define clade groups
df <- df %>%
  mutate(Clade_Group= case_when(
    Clade %in% c(5) ~ "Clade_5",
    Clade %in% c(1,2,3,4) ~ "Non_Clade_5",
    TRUE ~ NA_character_
  ))

# ignore control wells
df = df[df$Substrate %in% substrates,]
df = df[df$Clade_Group %in% c("Clade_5","Non_Clade_5"),]

# reduce to medians of technical replicates
df <- df %>%
  group_by(Substrate, Clade_Group, Substrate_Group, Isolate) %>%
  summarize(Carrying_Capacity = median(Carrying_Capacity, na.rm = TRUE), .groups = "drop")

# ensure that certain variables are treated as factors
df$Substrate <- as.factor(df$Substrate)
df$Substrate_Group <- as.factor(df$Substrate_Group)
df$Clade_Group <- as.factor(df$Clade_Group)

df$interaction <- (df$Substrate_Group == "Simple_Sugar") * (df$Clade_Group == "Clade_5")

# define linear mixed effects models
formula_full = 'Carrying_Capacity ~ interaction + Substrate_Group + Clade_Group + (1|Substrate) + (1|Isolate)'
formula_no_intxn = 'Carrying_Capacity ~ Substrate_Group + Clade_Group + (1|Substrate)  + (1|Isolate)'
formula_no_clade = 'Carrying_Capacity ~ Substrate_Group + (1|Substrate) + (1|Isolate)'

# run models
model_full = lmer(formula_full, data=df, REML=FALSE)
model_no_intxn = lmer(formula_no_intxn, data=df, REML=FALSE)
model_no_clade = lmer(formula_no_clade, data=df, REML=FALSE)
res_anova_intxn <- anova(model_full,model_no_intxn) %>% as.data.frame()
res_anova_clade <- anova(model_no_intxn,model_no_clade) %>% as.data.frame()

# display results
cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nFull Model')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(model_full)

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nRemove interaction term')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(model_no_intxn)
print(res_anova_intxn)

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nRemove clade group term')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(model_no_clade)
print(res_anova_clade)
cat('\n\n')


# store anova/tukey results in excel worksheet
worksheet <- "../../tables/linear_models_clade_5.xlsx"
write.xlsx(res_anova_clade, file=worksheet, sheetName="LMM-ANOVA-Clade", row.names=FALSE)
write.xlsx(res_anova_intxn, file=worksheet, sheetName="LMM-ANOVA-Interaction", row.names=FALSE,append=TRUE)

# question: does clade 5 grow differentially on each individual substrate?
# run two-sided Student's t-test

p_values <- numeric()
substrates <- unique(df$Substrate)
for (substrate in substrates) {

    df_stats = df[df$Substrate == substrate,]
    res <- df_stats %>% do(te=t.test(Carrying_Capacity ~ Clade_Group,alternative='two.sided',var.equal=TRUE,data=.))

    # get and store p-value
    p_value <- res$te[[1]]$p.value
    p_values <- c(p_values, p_value)
}

# FDR-correction of p-values
q_values <- p.adjust(p_values,method="fdr")

df_pvalues <- data.frame(
  Susbtrate = substrates, 
  p_value = p_values,
   q_value  = q_values
)

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nUnivariate test for each substrate')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(df_pvalues)
cat('\n\n')

write.xlsx(df_pvalues, file=worksheet, sheetName='T-tests-All', row.names=FALSE,append=TRUE)