#!/usr/bin/env Rscript

library(lme4)
library(dplyr)
library(xlsx)

# read growth data
df <- paste('../../amiga-ribotype-255/summary/merged_summary.txt')
df <- read.table(df,sep="\t",header=TRUE)
df <- df %>% filter(Isolate != "None" & Isolate != "")

# to harmonize analyses, rename column headers
names(df)[names(df) == 'Carbohydrate'] <- 'Substrate'
names(df)[names(df) == 'k_lin'] <- 'Carrying_Capacity'

# define substrates to analyse
substrates <- c("Fructose", "Ribose","None")
df$Substrate <- as.character(df$Substrate)
df <- df %>% filter(Substrate %in% substrates)

# define ribotype groups
df$Ribotype_Group = "Other"
df <- df %>% mutate(Ribotype_Group= ifelse(Ribotype == "RT255", "RT255", "Other"))

# reduce to medians of technical replicates
df <- df %>%
  group_by(Media, Substrate, Concentration_mM, Ribotype_Group, Isolate) %>%
  summarize(Carrying_Capacity = median(Carrying_Capacity, na.rm = TRUE), .groups = "drop")

# ensure that certain variables are treated as factors
df$Substrate <- as.factor(df$Substrate)
df$Ribotype_Group <- as.factor(df$Ribotype_Group)

# store anova results in excel worksheet
worksheet <- "../../tables/linear_models_ribotype_255.xlsx"

# run models
p_values <- numeric()
substrates <- c("Fructose", "Ribose")
for (substrate in substrates) {

    df_substrate = df[df$Media=="CDMM" & df$Substrate == substrate,]

    formula_full <- 'Carrying_Capacity ~ Ribotype_Group + Concentration_mM + (1|Isolate)'
    formula_null <- 'Carrying_Capacity ~ Concentration_mM + (1|Isolate)' 
    model_full = lmer(formula_full, data=df_substrate,REML=FALSE)
    model_null = lmer(formula_null, data=df_substrate,REML=FALSE)  
    res_anova <- anova(model_full,model_null) %>% as.data.frame() 

    cat('\n\n')
    print(substrate)
    cat('\n')
    print(model_full) 
    cat('\n')
    print(res_anova)

    # get and store p-value
    p_value <- res_anova$`Pr(>Chisq)`[2]
    p_values <- c(p_values, p_value)
    write.xlsx(res_anova, file=worksheet, sheetName=paste("LMM-ANOVA",substrate,sep="-"), row.names=FALSE,append=TRUE)
}

# FDR-correction of p-values
q_values <- p.adjust(p_values,method="fdr")

df_pvalues <- data.frame(
  Susbtrate = substrates, 
  p_value = p_values,
  q_value  = q_values
)

cat('\n\n')
print(df_pvalues)
cat('\n\n')

# question: does ribotype 255 grow to higher carrying capacity in each comparison?  
# run one-sided Student's t-test

p_values <- numeric()
s_values <- c()
c_values <- numeric()
substrates <- c("Fructose", "Ribose","None","BHIS")
for (substrate in substrates) {

    if (substrate != "BHIS") {
        df_substrate = df[df$Media=="CDMM" & df$Substrate == substrate,]
    } 
    else {
        df_substrate = df[df$Media=="BHIS",]
    }

    for (conc in unique(df_substrate$Concentration_mM)) {

        df_stats = df_substrate[df_substrate$Concentration_mM == conc,]
        res <- df_stats %>% do(te=wilcox.test(Carrying_Capacity ~ Ribotype_Group,alternative='less',correct=TRUE,data=.))
        
        # get and store p-value
        p_value <- res$te[[1]]$p.value
        p_values <- c(p_values, p_value)
        s_values = c(s_values,substrate)
        c_values = c(c_values,conc)
    }
}

# FDR-correction of p-values
q_values <- p.adjust(p_values,method="fdr")

df_pvalues <- data.frame(
  Susbtrate = s_values, 
  Concentration = c_values,
  p_value = p_values,
 q_value  = q_values
)

cat('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
cat("\nUnivariate test for each substrate")
cat('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
print(df_pvalues)
cat('\n\n')

write.xlsx(df_pvalues, file=worksheet, sheetName='Wilcoxon-All', row.names=FALSE,append=TRUE)