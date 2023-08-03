library(lme4)
library(rstatix)
library(rsq)
library(dplyr)

# Compare alpha for different conditons (does alpha change w/ set size or for guided compared to unguided?)
# Linear model 
dir_behav = "Z:\\Visual Search RFT\\results\\meg\\6 Alpha"

## Pre-search #####################################

file_name = "iaf_pow_bsl_condi.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$set.size <- factor(data$set.size)
data$guided <- factor(data$guided)


bsl_model <- lmer(iafpow~ 1 +(1|id),data=data)

model_set <- lmer(iafpow~set.size + (1|id),data=data)

model_set_gui <- lmer(iafpow~set.size + guided + (1|id),data=data)

anova(bsl_model,model_set,model_set_gui)

### Search ##############################################

file_name = "iaf_pow_search_condi.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$set.size <- factor(data$set.size)
data$guided <- factor(data$guided)
data$tp <- factor(data$tp)

bsl_model <- lmer(iafpow~ 1 +(1|id),data=data)

model_set <- lmer(iafpow~set.size + (1|id),data=data)

model_set_gui <- lmer(iafpow~set.size + guided + (1|id),data=data)

model_set_gui_tp <- lmer(iafpow~set.size + guided + tp + (1|id),data=data)


anova(bsl_model,model_set,model_set_gui,model_set_gui_tp)



