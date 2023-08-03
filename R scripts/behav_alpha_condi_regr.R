library(lme4)
library(rstatix)
library(rsq)
library(dplyr)
library(optimx)
library(reshape2)
#install.packages("optimx")



# Behavior for alpha high low (linear model)



dir_behav = "Z:\\Visual Search RFT\\results\\behavior"




## AVG RT Wilcox! #######################################

# change 250_500 to 1000_0 for alpha during search
file_name = "RT_avg_alpha250_500_ta_tp.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)


se_diff <- sd(data$RT_high - data$RT_low)/sqrt(31)



data <- melt(data)



data$variable <- factor(data$variable)

pwc = wilcox_test(data,value~variable,paired=TRUE,alternative = "less")
pwc

# Z-score
Z = qnorm(pwc$p/2)
Z


eff = wilcox_effsize(data,value~variable,paired=TRUE)
eff

Z = eff$effsize*sqrt(31)

## AVG hit rate Wilcox ######################

# change 250_500 to 1000_0 for alpha during search

file_name = "hit_avg_alpha250_500_ta_tp.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

se_diff <- sd(data$hit_high - data$hit_low)/sqrt(31)


data <- melt(data)


data$variable <- factor(data$variable)

pwc = wilcox_test(data,value~variable,paired=TRUE,alternative = "greater")
pwc

# Z-score



eff = wilcox_effsize(data,value~variable,paired=TRUE)
eff

Z = eff$effsize*sqrt(31)




## RT LME ###############################################


file_name = "RT_hit_alpha250_500_condi_ta_tp.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$id <- factor(data$id)

# decode condition into set size and guided/unguided

data$gui.ung[data$condi==1] <- 0
data$gui.ung[data$condi==2] <- 1
data$gui.ung[data$condi==3] <- 0
data$gui.ung[data$condi==4] <- 1

data$gui.ung[data$condi==5] <- 0
data$gui.ung[data$condi==6] <- 1
data$gui.ung[data$condi==7] <- 0
data$gui.ung[data$condi==8] <- 1


data$set.size[data$condi==1] <- 1
data$set.size[data$condi==3] <- 2
data$set.size[data$condi==2] <- 1
data$set.size[data$condi==4] <- 2

data$set.size[data$condi==5] <- 1
data$set.size[data$condi==7] <- 2
data$set.size[data$condi==6] <- 1
data$set.size[data$condi==8] <- 2


data$tp[data$condi %in% c(1,2,3,4)] <- 0
data$tp[data$condi %in% c(5,6,7,8)] <- 1



## RT #################################


RT_high = data$RT[data$high_low==1]
RT_low = data$RT[data$high_low==0]

RT_diff = RT_high - RT_low

shapiro.test(RT_diff)


hist(RT_diff)
dens <- dnorm(min(RT_diff),max(RT_diff))
plot(density(RT_diff))




# models
bsl_model <- glmer(RT~ 1 +(1|id),data=data, family=Gamma,control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                                                optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

model_setsize <- glmer(RT~1 + set.size + (1 +set.size|id),data=data, family=Gamma(identity),
                       control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                              optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))


model_setsize_gui <- glmer(RT~ 1+ set.size + gui.ung + (1 +set.size|id) + (1 +gui.ung|id),data=data, family=Gamma(identity),
                           control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                  optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

model_setsize_gui_tp <- glmer(RT~set.size + gui.ung + tp + (1 +set.size|id) + (1 +gui.ung|id) + (1 +tp|id),data=data, family=Gamma(identity),
                              control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                     optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))

model_setsize_gui_tp_alpha <- glmer(RT~set.size + gui.ung + tp + high_low + (1 +set.size|id) + (1 +gui.ung|id) + (1 +tp|id) + (1 +high_low|id),data=data, family=Gamma(identity),
                                    control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                           optCtrl = list(method = "L-BFGS-B", starttests = FALSE, kkt = FALSE)))




anova(bsl_model,model_setsize,model_setsize_gui,model_setsize_gui_tp,model_setsize_gui_tp_alpha)

rsq(bsl_model)

rsq(model_setsize)
rsq(model_setsize_gui)

rsq(model_setsize_gui_tp)
rsq(model_setsize_gui_tp_alpha)

rsq(model_setsize_gui_tp_alpha)$model - rsq(model_setsize_gui_tp)$model


summary(model_setsize_gui_tp_alpha)


# DF for fixed effects
n <- nrow(model_setsize_gui@frame)
p <- length(coef(summary(model_setsize_gui_tp_alpha)))
df <- n-p


# Wilcoxon
pwc <- data %>%
  group_by(set.size,gui.ung,tp) %>%
  pairwise_wilcox_test(RT ~ high_low,paired=TRUE,alternative = "greater") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
pwc


# SE of difference

# standard error of the difference
data_high = data[data$high_low==1,]
data_low = data[data$high_low==0,]

RT_highvslow = data_high
RT_highvslow$diff = data_high$RT - data_low$RT

rt_diff_sum <- RT_highvslow %>%
  group_by(set.size,gui.ung,tp) %>%
  summarise(
    Weight = mean(diff),
    sd = sd(diff),
    n = n(),
    se = sd / sqrt(n)
  )


pwc <- data %>%
  group_by(set.size,gui.ung,tp) %>%
  wilcox_effsize(RT ~ high_low,paired=TRUE,alternative = "greater")

pwc




## Hit rate #################################

# models
bsl_model <- glmer(d~1 +(1|id),data=data)
model_setsize <- glmer(d~1 + set.size + (1 +set.size|id),data=data)
model_setsize_gui <- glmer(d~1+ set.size + gui.ung + (1 +set.size|id) + (1 +gui.ung|id),data=data)

model_setsize_gui_tp <- glmer(hit.rate~set.size + gui.ung + tp + + (1 +set.size|id) + (1 +gui.ung|id) + (1 +tp|id),data=data)
model_setsize_gui_tp_alpha <- lmer(hit.rate~set.size + gui.ung + tp + + high_low + (1 +set.size|id) + (1 +gui.ung|id) + (1 +tp|id) + (1 +high_low|id),data=data)



anova(bsl_model,model_setsize,model_setsize_gui,model_setsize_gui_tp,model_setsize_gui_tp_alpha)

rsq(bsl_model)

rsq(model_setsize)
rsq(model_setsize_gui)

rsq(model_setsize_gui_tp)
rsq(model_setsize_gui_tp_alpha)


summary(model_setsize_gui_tp_alpha)
