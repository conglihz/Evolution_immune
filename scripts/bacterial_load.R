setwd("~/Documents/Zhao_Lab/manuscript/figures/bacterial_load")
library(readxl)
library(ggplot2)
library(ggsci)
library(cowplot)

data <- read_excel("bacterial_load_1day.xlsx")

head(data)

data$Species <- factor(data$Species, levels = c('dmel', 'dsim', 'dana', 'dvir', 'sleb'))

data_summary <- data %>%
  group_by(Species) %>%
  summarise(
    median_load = median(Load, na.rm = TRUE),
    mean_load = mean(Load, na.rm = TRUE),
    se_load = sd(Load, na.rm = TRUE) / sqrt(n())
  )


library(dplyr)
library(tidyr)
library(purrr)

# Function to return 95% CI and median
bootstrap_median_ci <- function(x, n_boot = 1000, conf = 0.95) {
  medians <- replicate(n_boot, median(sample(x, replace = TRUE), na.rm = TRUE))
  quantile(medians, probs = c((1 - conf)/2, 0.5, 1 - (1 - conf)/2), na.rm = TRUE, names = FALSE)
}
data_summary <- data %>%
  group_by(Species) %>%
  summarise(ci = list(bootstrap_median_ci(Load))) %>%
  unnest_wider(ci, names_sep = "_") %>%
  rename(
    ci_lower = ci_1,
    median_load = ci_2,
    ci_upper = ci_3
  )


# find default log10 breaks
default_breaks <- scales::log_breaks()(range(data$Load, na.rm = TRUE))

# add extra tick at 0.0001
custom_breaks <- sort(unique(c(0.01, default_breaks)))

p <- ggplot(data, aes(x = Species, y = Load, color = Species)) +
      # plot individual points (jittered)
      geom_point(
        size = 0.3,
        alpha = 1,  # lower alpha for individual points
        position = position_jitter(width = 0.2)  # add jitter to avoid overlap
      ) +
      # vertical CI bars
      geom_errorbar(
        data = data_summary,
        aes(x = Species, ymin = ci_lower, ymax = ci_upper),
        width = 0.2,
        size = 0.3,
        color = "black",
        inherit.aes = FALSE
      ) +
      
      # horizontal line at median
      geom_errorbar(
        data = data_summary,
        aes(x = Species, ymin = median_load, ymax = median_load),
        width = 0.5,
        size = 0.3,
        color = "black",
        inherit.aes = FALSE
      ) + 
      labs(
        title = "Bacterial load on day1",
        x = NULL,
        y = "CFU/ml",
        color = "Species"
      ) +
      scale_y_continuous(
        trans = "log10",
        breaks = custom_breaks,
        labels = function(x) {
          x_lab <- as.character(x)
          x_lab[x == 0.0001] <- "0"
          x_lab
        }
      ) +
      scale_color_npg() +  # apply NPG color palette
      theme(
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7, color = "black"),
        legend.position = "none",
        axis.line = element_line(color = "black", size = 0.4),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank()
      )
print(p)
p <- plot_grid(p, labels = 'C')
ggsave("bacterial_load.pdf", p, width = 2, height = 2, dpi = 300)

