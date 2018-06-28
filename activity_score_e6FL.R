setwd('C:/Users/lkoneva.UMROOT/Documents/HNSCC/E6E7_Yidan/lkoneva/both_results/activity_socre/e6ratio')
dir()

library(dplyr)
require(gridExtra)
library(matrixStats) # for score calculation on matrix
library(ggplot2)

df <- read.delim("both_voom_e6FL_16_gene_list_nooutl_model3cov.txt")
head(df)

sign <- df %>% filter(FDR_e6FL < 0.05)
head(sign)
#write.csv(sign, "both_voom_e6FL_16_nooutl_model3cov_sign.csv")

sign_pos <- df %>% filter(FDR_e6FL < 0.05 & slope_e6FL > 0)
sign_neg <- df %>% filter(FDR_e6FL < 0.05 & slope_e6FL < 0)

names(sign_pos)
genes_pos <- sign_pos[, 14:81]
row.names(genes_pos) <- sign_pos$gene_symbol
tail(genes_pos)

# transpose matrix
x <- as.data.frame(t(genes_pos))
head(x)

# rank samples by expression in ascending order for positively assotiated genes (coefficient > 0)
rank_pos <- as.matrix(mutate_all(x, funs(min_rank)))
head(rank_pos)
row.names(rank_pos) <- row.names(x)
head(rank_pos[, 1:6],20) # check correct order of ranks

# ranking for negatively associated genes
head(sign_neg)
genes_neg <- sign_neg[, 14:81]
row.names(genes_neg) <- sign_neg$gene_symbol
y <- as.data.frame(t(genes_neg))
head(y)

# rank samples by expression in descending order for negatively assotiated genes (coefficient < 0)
#rank_neg <- as.matrix(mutate_all(y, funs(dense_rank(desc(.)))))
rank_neg <- as.matrix(mutate_all(y, funs(min_rank(desc(.)))))
row.names(rank_neg) <- row.names(y)
head(rank_neg)
head(rank_neg[, 1:6],20) # check correct order of ranks

# merge rank_pos and rank_neg
rank_all <- merge(rank_pos, rank_neg, by = "row.names")
tail(rank_all[,1:6],20)

# convert to matrix
rank_mt <- as.matrix(rank_all[,-1])
head(rank_mt[,1:6])
row.names(rank_mt) <- rank_all$Row.names

# centering and scalin
rank_score <- transform(rank_mt, sum = rowSums(rank_mt), mean = rowMeans(rank_mt), SD = rowSds(rank_mt))

# calculating scores
rank_score2 <- as.data.frame(rank_score) %>% mutate(score = (sum-mean)/SD)
rank_score2$sample <- row.names(rank_mt)
head(rank_score2)
#write.csv(rank_score2, "activity_score_e6FL.csv", row.names = F)

##################################################

# read meta file prepared for correlation clinical with activiti scores
meta <- read.csv("meta_both_cohorts_hpv16_for_activity_score.csv")

# remove outlier sample TCGA-CN-A6V1 from meta
meta
meta.nooutl <- meta %>% filter(ID != 'TCGA-CN-A6V1')

# select columns (ID_3, Sex, Age, Smoking.Status, HPV.integration, TUMOR_CLASSIF_T_clean, Cohort,
# Site_combned, Stage_combined, Nodal01vs23, e6_fl_ratio)
names(meta.nooutl)
meta2 <- meta.nooutl[,c(3,4,5,6,7,9,10,11,12,15,16,17,18,19,32)]
head(meta2)

# Continue here: bind with meta file
head(rank_score2)
score <- rank_score2[, c(174,173)]
tail(score)
meta_score <- merge(score, meta2, by.x="sample", by.y="ID_3")
head(meta_score)

#write.csv(meta_score, "meta_score_e6FL.csv")

# perform anova on scores between clinical groups: T-stage etc.
#meta_score <- read.csv("meta_score_e6FL.csv")
names(meta_score)

# score vs e6_fl_ratio - to check that score correlates with efFL
fit = lm(score ~ e6_fl_ratio, meta_score)
fit
anova(fit)
signif(summary(fit)$adj.r.squared, 3)
signif(summary(fit)$coef[2,4], 3)

results <- data.frame()
for(i in 3:16){
  fit = lm(score ~ meta_score[,i], meta_score)
#  r.adj <- signif(summary(fit)$adj.r.squared, 3)
#  pval <- signif(summary(fit)$coef[2,4], 3)
  newline <- data.frame(Parameter = colnames(meta_score[i]), 
                        Adj.R2 = signif(summary(fit)$adj.r.squared, 3), 
                        P.value = signif(summary(fit)$coef[length(summary(fit)$coefficients[,4]),4], 3),
                        Intercept = signif(fit$coef[[1]],3 ),
                        Slope = signif(fit$coef[[2]], 3))
  results <- rbind(results, newline) 
}
results
write.csv(results, "correlation_activity_score_E6FL_clinical.csv", row.names = F)

##### directions
cor.test(meta_score$score, meta_score$e6_fl_ratio)

table(meta_score$Site_combined)
summary(meta_score$score[meta_score$Site_combined == "Oropharynx"])
summary(meta_score$score[meta_score$Site_combined == "Other"])

table(meta_score$TUMOR_CLASSIF_T_clean)
summary(meta_score$score[meta_score$TUMOR_CLASSIF_T_clean == "1"])
summary(meta_score$score[meta_score$TUMOR_CLASSIF_T_clean == "2"])
summary(meta_score$score[meta_score$TUMOR_CLASSIF_T_clean == "3"])
summary(meta_score$score[meta_score$TUMOR_CLASSIF_T_clean == "4"])

table(meta_score$HPV.integration)
summary(meta_score$score[meta_score$HPV.integration == "neg"])
summary(meta_score$score[meta_score$HPV.integration == "pos"])

table(meta_score$Smoking.Status)
summary(meta_score$score[meta_score$Smoking.Status == "Current"])
summary(meta_score$score[meta_score$Smoking.Status == "Former"])
summary(meta_score$score[meta_score$Smoking.Status == "Never"])

########
meta_score
meta_score$T.stage <- as.factor(meta_score$TUMOR_CLASSIF_T_clean)

png("boxplot_score_T-stage_withNA.png",width = 8, height = 6, units = 'in', res = 300)
p = ggplot(meta_score, aes(x=T.stage, y=score, fill=T.stage)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(title="Activity score by tumor's T-stages",x="T-stages", y = "Activity score")
p + theme_classic()
dev.off()

# remove one sample with T-stage = NA
meta_score2 <- meta_score %>% filter(is.na(TUMOR_CLASSIF_T_clean) == FALSE)

png("boxplot_score_T-stage.png",width = 8, height = 6, units = 'in', res = 300)
p = ggplot(meta_score2, aes(x=T.stage, y=score, fill=T.stage)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(title="Activity score by tumor's T-stages",x="T-stages", y = "Activity score")
p + theme_classic()
dev.off()

head(meta_score)

png("boxplot_scoree6ratio_Site_combined.png",width = 8, height = 6, units = 'in', res = 300)
p = ggplot(meta_score, aes(x=Site_combined, y=score, fill=Site_combined)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(title="Activity score by tumor's Site (combined)",x="Tumor site", y = "Activity score")
p + theme_classic()
dev.off()

png("boxplot_scoree6ratio_HPVintegration.png",width = 8, height = 6, units = 'in', res = 300)
p = ggplot(meta_score, aes(x=HPV.integration, y=score, fill=HPV.integration)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(title="Activity score by tumor's HPV-integration",x="HPV integration", y = "Activity score")
p + theme_classic()
dev.off()

png("boxplot_scoree6ratio_Smoking_Status.png",width = 8, height = 6, units = 'in', res = 300)
p = ggplot(meta_score, aes(x=Smoking.Status, y=score, fill=Smoking.Status)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.1)) +
  labs(title="Activity score by smoking status",x="Smoking status", y = "Activity score")
p + theme_classic()
dev.off()

########
# look at ANOVA statistics by separate variables
# Sex Age Site Stage_clean Smoking.Status HPV.integration TUMOR_CLASSIF_T_clean Nodal.Status_clean cohort
# Site_combined Stage_combined Nodal_combined Nodal01vs23 loge6e7

# score vs var
var = meta_score$Nodal01vs23
table(var)
fit = lm(score ~ var, meta_score)
fit
anova(fit)
summary(fit)

# consider T-stage as factor with 4 levels not as numerical variable
var = as.factor(meta_score$TUMOR_CLASSIF_T_clean)
table(var)
str(var)
fit = lm(score ~ var, meta_score)
fit
anova(fit)
summary(fit)

# exclude NA
meta_score3 <- meta_score[-53, ]
var3 = meta_score3$Nodal.Status_clean
table(var3)
str(var3)
fit3 = lm(score ~ var3, meta_score3)
fit3
anova(fit3)
summary(fit3)
