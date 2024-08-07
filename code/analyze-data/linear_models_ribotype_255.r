#!/usr/bin/env Rscript

library(lme4)
library(dplyr)
library(xlsx)

# read growth data
df <- paste('../../amiga-ribotype-255/summary/merged_summary.txt')
df <- read.table(df,sep="\t",header=TRUE)

# to harmonize analyses, rename column headers
names(df)[names(df) == 'Carbohydrate'] <- 'Substrate'
names(df)[names(df) == 'k_lin'] <- 'Carrying_Capacity'

# define substrates to analyse
substrates <- c("Fructose", "Ribose","None")

# define ribotype groups
df$Ribotype_Group = "Other
"
df <- df %>%
  mutate(Ribotype_Group= ifelse(Ribotype == "RT255", "RT255", "Other")
  )

# question: does ribotype 255 grow differentially regardless of media and substrate type?

formula_full = 'Carrying_Capacity ~ Ribotype_Group*Substrate + Media + (1|Isolate)'
formula_no_intxn = 'Carrying_Capacity ~ Ribotype_Group + Substrate + Media + (1|Isolate)'
formula_no_clade = 'Carrying_Capacity ~ Substrate + Media + (1|Isolate)'

model_full = lmer(formula_full, data=df, REML=FALSE)
model_intxn = lmer(formula_no_intxn, data=df, REML=FALSE)
model_clade = lmer(formula_no_clade, data=df, REML=FALSE)
res_anova_intxn <- anova(model_full,model_intxn) %>% as.data.frame()
res_anova_clade <- anova(model_full,model_clade) %>% as.data.frame()

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nFull Model')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(model_full)

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nRemove interaction term')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(model_intxn)

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat('\nRemove ribotype group term')
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(model_clade)
cat('\n\n')

# store anova results in excel worksheet
worksheet <- "../../tables/linear_models_ribotype_255.xlsx"
write.xlsx(res_anova_intxn, file=worksheet, sheetName="main_intxn", row.names=FALSE)
write.xlsx(res_anova_clade, file=worksheet, sheetName="main_clade", row.names=FALSE,append=TRUE)

# question: does ribotype 255 grow differentially on each individual substrate?

formula_full = 'Carrying_Capacity ~ Ribotype + (1|Isolate)'
formula_null = 'Carrying_Capacity ~ (1|Isolate)'

p_values <- numeric()

for (substrate in substrates) {

    df_substrate = df[df$Media=="CDMM",]
    df_substrate = df_substrate[df_substrate$Substrate==substrate,]

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