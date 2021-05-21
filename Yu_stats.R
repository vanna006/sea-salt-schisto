#####################################
# Stastistics for Yu et al. Sea salt#
#      Parasitism interactions      #
#####################################

setwd('E:/PhD/Undergrads/Anna_Ao_Yu')

# examine snail reproductive output
eggs <- read.csv('egg.csv')
hist(eggs$eggs, breaks = 15)

library(ggpubr)
library(glmmTMB)
library(bbmle)
library(car)
library(emmeans)

#a subset of only control snails
green <- eggs[1:480,]
green.no <- na.omit(green)

#for green infections
ggboxplot(green.no, x = 'week', y = 'eggs', color = 'fact')
fit_zipoisson <- glmmTMB(eggs ~ fact * week + (1|ID), data=green.no, ziformula=~1, family=poisson)
fit_hzip <- update(fit_zipoisson,
                   ziformula = ~.,
                   data = green.no,
                   family = poisson(link="log"))
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_hzip, family=list(family="truncated_nbinom1",link="log"))
AICtab(fit_zipoisson,fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)
summary(fit_hzip)
Anova(fit_hzip, type='III')

# green reproduction was different with infection so need to look specifically
# at the uninfected snails or infected

uninf <- subset(eggs, infection < 1)
ggboxplot(uninf, x = 'week', y = 'eggs', color = 'group')

#generate a general model and create models of various distribution
fit_zipoisson <- glmmTMB(eggs ~ group + (1|ID), data=uninf, ziformula=~1, family=poisson)
fit_hzip <- update(fit_zipoisson,
                   ziformula = ~.,
                   data = uninf,
                   family = poisson(link="log"))
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_hzip, family=list(family="truncated_nbinom1",link="log"))

#determine the best fitting model
AICtab(fit_zipoisson,fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)

#model output and statistics
summary(fit_zinbinom1)
Anova(fit_zinbinom1, type='III')
multzip <- emmeans(fit_zinbinom1, pairwise ~ group, adjust = 'mvt')
multzip

#generate violin plot of data
uninf.g <- ggplot(uninf, aes(x = gfct, y = eggs, fill = group)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Egg masses", hjust=10) +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle=90, hjust=1))

saltlabels <- c('0.38 ppt', '1.9 ppt', '3.8 ppt', '5.7 ppt')

uninf.g + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 0.5, position=position_dodge(1)) +
  scale_fill_manual(values=c("green4","orangered" ,'blue4' , 'yellow3')) + theme_classic() +
  theme(text = element_text(size=10), legend.position='none') + scale_x_discrete(labels= saltlabels)

######### Examine reproduciton amongst infected classes
infec <- subset(eggs, infection > 0)
ggboxplot(infec, x = 'group', y = 'eggs', color = 'group')
fit_zipoisson <- glmmTMB(eggs ~ group + (1|ID), data=infec, ziformula=~1, family=poisson)
fit_hzip <- update(fit_zipoisson,
                   ziformula = ~.,
                   data = infec,
                   family = poisson(link="log"))
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_hzip, family=list(family="truncated_nbinom1",link="log"))
AICtab(fit_zipoisson,fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)
summary(fit_hnbinom2)
Anova(fit_hnbinom2, type='III')
multzip <- emmeans(fit_hnbinom2, pairwise ~ group, adjust = 'mvt')
multzip

###### Snail survival analysis
library("survival")
library("survminer")
salty <- read.csv("anna.surv.csv")
fit <- survfit(Surv(time, status) ~ group, data = salty)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = F, conf.int = TRUE,
           palette = c("green4", "orangered", 'blue4', 'yellow3'),
           xlab = "Week",
           #risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

surv_diff <- survdiff(Surv(time, status) ~ group, data = sun)
surv_diff

cph <- coxph(Surv(time, status) ~ group, data = salty)
cph
summary(cph)
multzip <- emmeans(cph, list(pairwise ~ group), adjust = 'mvt')
multzip


psurv_diff <- pairwise_survdiff(Surv(time, status) ~ group, data = salty, p.adjust.method = 'hochberg')
psurv_diff


## Look fo differences in survival between infected and uninfected over experiment
gr <- salty[1:60,]
gr.no <- na.omit(gr)

fitgr <- survfit(Surv(time, status) ~ fact, data = gr.no)
print(fitgr)
summary(fitgr)

ggsurvplot(fitgr,
           pval = F, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

surv_diff <- survdiff(Surv(time, status) ~ fact, data = gr.no)
surv_diff

cph <- coxph(Surv(time, status) ~ fact, data = gr.no)
cph
summary(cph)

psurv_diff <- pairwise_survdiff(Surv(time, status) ~ fact, data = gr.no, p.adjust.method = 'hochberg')
psurv_diff



### Examine differences in prevalence across groups
salty.no <- na.omit(salty)
prev <- glm(infection ~ group, data = salty, family = 'binomial')
summary(prev)
multzip <- emmeans(prev, pairwise ~ group, adjust = 'mvt')
multzip

H<- c(0.702, 0.161, 0.132, 0.093)

barplot(H)
par(mar=c(8.1, 5, 4.1, 2.1))
barplot(H, xlab = " 0.38 ppt       1.9 ppt         3.8 ppt         5.7 ppt", ylab = "Proportion of hosts infected", ylim = c(0,1), cex.axis = 1, cex.lab= 1,
        col = c("green4", 'yellow3', "orangered", 'blue4'))
abline(h=0)

### Cercarial survival analysis; have to do mixed effects survival analysis with coxme
library(coxme)
cercs.a <- read.csv("csurv.ao.csv")

#generate a standard model in order to plot data
coxfit <- survfit(Surv(time, status) ~ confct, data = cercs.a)

#generate mixed effects cox model for stats
cox.m <- coxme(Surv(time, status) ~ confct + (1|snailID), data = cercs.a)
cox.m

ggsurvplot(coxfit,
           pval = F, conf.int = TRUE,
           xlab = "Week",
           palette = c("green4", 'yellow3', "orangered", 'blue4'),
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

multzip <- emmeans(cox.m, pairwise ~ confct, adjust = 'mvt')
multzip










