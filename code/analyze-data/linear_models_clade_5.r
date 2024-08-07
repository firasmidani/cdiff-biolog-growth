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

# question: does substrate group differentially impact growth by clade 5?

formula_full = 'Carrying_Capacity ~ Clade_Group*Substrate_Group + Substrate + (1|Isolate)'
formula_no_intxn = 'Carrying_Capacity ~ Clade_Group + Substrate_Group + Substrate + (1|Isolate)'
formula_no_clade = 'Carrying_Capacity ~ Substrate_Group + Substrate + (1|Isolate)'

model_full = lmer(formula_full, data=df, REML=FALSE)
model_no_intxn = lmer(formula_no_intxn, data=df, REML=FALSE)
model_no_clade = lmer(formula_no_clade, data=df, REML=FALSE)
res_anova_intxn <- anova(model_full,model_no_intxn) %>% as.data.frame()
res_anova_clade <- anova(model_no_intxn,model_no_clade) %>% as.data.frame()

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nFull Model')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(model_full)

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nRemove interaction term')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(model_no_intxn)

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nRemove clade group term')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(model_no_clade)
cat('\n\n')

# store anova results in excel worksheet
worksheet <- "../../tables/linear_models_clade_5.xlsx"
write.xlsx(res_anova_intxn, file=worksheet, sheetName="main_intxn", row.names=FALSE)
write.xlsx(res_anova_clade, file=worksheet, sheetName="main_clade", row.names=FALSE,append=TRUE)

# question: does clade 5 grow differentially on each individual substrate?

formula_full = 'Carrying_Capacity ~ Clade_Group + (1|Isolate)'
formula_null = 'Carrying_Capacity ~ (1|Isolate)'

p_values <- numeric()

for (substrate in substrates) {

    df_substrate = df[df$Substrate==substrate,]

    model_full = lmer(formula_full, data=df_substrate, REML=FALSE)
    model_null = lmer(formula_null, data=df_substrate, REML=FALSE)
    res_anova <- anova(model_full,model_null) %>% as.data.frame()

    # get and store p-value
    p_value <- res_anova$`Pr(>Chisq)`[2]
    p_values <- c(p_values, p_value)

    write.xlsx(res_anova, file=worksheet, sheetName=substrate, row.names=FALSE,append=TRUE)
}

# FDR-correction of p-values
q_values <- p.adjust(p_values,method="fdr")

df_pvalues <- data.frame(
  Susbtrate = substrates, 
  p_value = p_values,
   q_value  = q_values
)

print(df_pvalues)
cat('\n\n')