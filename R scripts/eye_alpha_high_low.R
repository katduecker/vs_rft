library(lme4)
library(rstatix)
library(rsq)
library(dplyr)
dir_behav = "Z:\\Visual Search RFT\\results\\eyelink"

# Compare ocular artefacts for alpha high vs low

## Blinks ##########################


# change 250_500 to 1000_0 for alpha during search
file_name = "blinks_alpha_high_low_250_500.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$id <- factor(data$id)

# decode condition into set size and guided/unguided

data$gui.ung[data$condi==1] <- 0
data$gui.ung[data$condi==2] <- 1
data$gui.ung[data$condi==3] <- 0
data$gui.ung[data$condi==4] <- 1

data$set.size[data$condi==1] <- 1
data$set.size[data$condi==3] <- 2
data$set.size[data$condi==2] <- 1
data$set.size[data$condi==4] <- 2

# models
bsl_model <- lmer(num_bl~ 1 +(1|id),data=data)
model_high<- lmer(num_bl~high + (1|id),data=data)
anova(bsl_model,model_high)

rsq(model_high)$model - rsq(bsl_model)$model 

## Saccades ##########################

file_name = "saccades_alpha_high_low.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$id <- factor(data$id)

# decode condition into set size and guided/unguided

data$gui.ung[data$condi==1] <- 0
data$gui.ung[data$condi==2] <- 1
data$gui.ung[data$condi==3] <- 0
data$gui.ung[data$condi==4] <- 1

data$set.size[data$condi==1] <- 1
data$set.size[data$condi==3] <- 2
data$set.size[data$condi==2] <- 1
data$set.size[data$condi==4] <- 2

# models
bsl_model <- lmer(num_sac~ 1 +(1|id),data=data)
model_high<- lmer(num_sac~high + (1|id),data=data)

anova(bsl_model,model_high)


rsq(model_high)$model - rsq(bsl_model)$model

# t-test
pwc <- data %>%
  group_by(gui.ung,set.size) %>%
  pairwise_t_test(num_sac ~ high, paired=TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
pwc

pwc <- data %>%
  group_by(gui.ung,set.size) %>%
  cohens_d(num_sac ~ high, paired=TRUE,hedges.correction = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

## Gaze bias #########################

file_name = "bias_alpha_high_low.csv"

file_path = file.path(dir_behav,file_name)
data <- read.csv(file_path)

data$id <- factor(data$id)

# decode condition into set size and guided/unguided

data$gui.ung[data$condi==1] <- 0
data$gui.ung[data$condi==2] <- 1
data$gui.ung[data$condi==3] <- 0
data$gui.ung[data$condi==4] <- 1

data$set.size[data$condi==1] <- 1
data$set.size[data$condi==3] <- 2
data$set.size[data$condi==2] <- 1
data$set.size[data$condi==4] <- 2


bsl_model <- lmer(bias~ 1 +(1|id),data=data)
model_high<- lmer(bias~high + (1|id),data=data)

anova(bsl_model,model_high)
rsq(model_high)$model - rsq(bsl_model)$model
