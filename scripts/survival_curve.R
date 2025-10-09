setwd("~/Documents/Zhao_Lab/infection/survival")
library(survival)
library(survminer)
library(dplyr)
library(ggsci)
library(cowplot)
library(ggplot2)
library(ggrepel)

# Pr
data = read.csv('pr_data.csv')
glimpse(data)

sps=c('dmel','dsim','dana','dvir','sleb')
for (s in sps)
{
  data2 = data[data$sp == s,]
  t = paste("Dead_day distribution in",s)
  hist(data2$dead_day, main = t, xlim= c(0,30), ylim = c(0,50), breaks = 30)
}

#transform grouping variable(sp in this case) to factor class so I can define levels
data1 <- data %>% mutate(sp= factor(sp, levels = c('dmel','dsim','dana','dvir','sleb')))

# fit survival data using the Kaplan-Meier method
surv_object <- Surv(time = data1$dead_day, event = data1$eventtype)
surv_object 

fit1 <- survfit(surv_object ~ sp, data = data1)
summary(fit1)

# get median survival
get_median_survival <- function(fit) {
  surv_summ <- summary(fit)
  result <- list()
  
  # get strata info
  strata_names <- names(fit$strata)
  strata_lengths <- fit$strata
  
  start <- 1
  for (i in seq_along(strata_names)) {
    end <- start + strata_lengths[i] - 1
    group_surv <- surv_summ$surv[start:end]
    group_time <- surv_summ$time[start:end]
    
    #remove NA values
    valid_idx <- which(!is.na(group_surv) & !is.na(group_time))
    
    if (length(valid_idx) == 0) {
      median_time <- NA
    } else {
      group_surv <- group_surv[valid_idx]
      group_time <- group_time[valid_idx]
      
      # find first time survival drops to or below 0.5
      if (any(group_surv <= 0.5)) {
        median_time <- group_time[min(which(group_surv <= 0.5))]
      } else {
        median_time <- NA  # median not reached
      }
    }
    
    result[[strata_names[i]]] <- median_time
    start <- end + 1
  }
  
  return(result)
}

get_median_survival(fit1)

# extract summary at day 3
summary_fit_day3 <- summary(fit1, times = 3)

# create data frame
survival_day3_df <- data.frame(
  species = gsub("sp=", "", summary_fit_day3$strata),
  time = summary_fit_day3$time,
  n_risk = summary_fit_day3$n.risk,
  n_event = summary_fit_day3$n.event,
  survival = summary_fit_day3$surv,
  std_err = summary_fit_day3$std.err,
  lower_ci = summary_fit_day3$lower,
  upper_ci = summary_fit_day3$upper
)
write.csv(survival_day3_df, file = "survival_day3.csv", row.names = FALSE)





# Fit a Cox proportional hazards model
fit.coxph <- coxph(surv_object ~ sp + batch, 
                   data = data1)
ggforest(fit.coxph, data = data1)

fit.coxph


#3day cox
data_3days <- subset(data1, dead_day <= 3)
surv_object_3day <- Surv(time = data_3days$dead_day, event = data_3days$eventtype)

fit.coxph <- coxph(surv_object_3day ~ sp + batch, 
                   data = data_3days)
ggforest(fit.coxph, data = data_3days)


summary(fit1, times = 7)

p <- ggsurvplot(fit1, data = data1, palette = 'npg',
           title = "Providencia rettgeri",
           pval = FALSE, conf.int = FALSE, 
           xlab = 'Time (days)', xlim = c(0,29), break.x.by = 5,
           ylab = 'Survival probability (%)',
           fun = "pct", # make y axis percentage
           legend = "right",
           legend.title = "Species",
           legend.labs= c('dmel','dsim','dana','dvir','sleb'),
           size=0.5
           )
####################
# style1 with legend
pr_plot <- p$plot +
  ggtitle(expression(italic("Providencia rettgeri"))) +
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        axis.title.x = element_text(size = 8),       # font.x
        axis.title.y = element_text(size = 8),       # font.y
        plot.caption = element_text(size = 8),       # font.caption
        legend.title = element_text(size = 7),       # font.legend (legend title)
        legend.text = element_text(size = 7),        # font.legend (legend labels)
        axis.text.x = element_text(size = 8),        # font.tickslab (x-axis tick labels)
        axis.text.y = element_text(size = 8),         # font.tickslab (y-axis tick labels)
        axis.line = element_line(size = 0.4),
        axis.ticks = element_line(size = 0.4)
        #text = element_text(size = rel(1))
  )
print(pr_plot)

ggsave("Pr_survival_0213.pdf", pr_plot, width = 3.5, height = 3, dpi = 300)

