library(ggplot2)
library(tidyverse)

here::i_am("capture_efficiency.R")

dpi <- 300

# Count per spot

df <- read.csv(here::here("count_per_spot.csv")) %>%
  mutate(condition = factor(condition, levels = c("9m_AD_R1", "9m_AD_R2", "9m_WT_R1", "9m_WT_R2"))) %>%
  filter(count > 0)

medians <- df %>%
  group_by(condition) %>%
  summarise(med = median(count)) %>%
  mutate(condition_label = paste0(condition, "\nmedian: ", round(med, 1))) %>%
  mutate(condition_label = factor(condition_label, levels = condition_label))

df_annotated <- df %>%
  left_join(medians, by = "condition")

p <- ggplot(df_annotated, aes(x = count)) +
  geom_histogram(binwidth = 500, color = "#e9ecef", fill = "#6fafd2", alpha = 0.9) +
  geom_vline(data = medians, aes(xintercept = med), color = "red", linetype = "dashed") +
  labs(x = "mRNA Count", y = "Grid Count") +
  facet_wrap(~ condition_label, nrow = 1, ncol = 4) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 15))

ggsave(here::here("count_per_spot.jpeg"), p, width = 16, height = 4, dpi = dpi)

# Gene capture efficiency

df <- read.csv(here::here("gene_capture_efficiency.csv")) %>%
  mutate(condition = factor(condition, levels = c("9m_AD_R1", "9m_AD_R2", "9m_WT_R1", "9m_WT_R2")))

ref_genes <- c('Bsn', 'Ank3', 'Syp', 'Slc32a1', 'Nav1', 'Map2', 'Gria1', 'Map1a', 'Homer1', 'Shank3', 'Nfasc', 'Gphn', 'Tubb3', 'Mapt', 'Slc17a6', 'Vamp2', 'Stx1a', 'Dlg3', 'Camk2a', 'Slc17a7', 'Gria2', 'Cplx2', 'Dlg4', 'Nlgn2', 'Gap43', 'Nrxn1', 'Homer2', 'Cyfip2', 'Syt1', 'Nlgn3', 'Ddn', 'Shank1', 'Nlgn1', 'Syn1')
df <- df[df$genes %in% ref_genes, ]

compare_groups <- c("9m_AD_R1", "9m_AD_R2", "9m_WT_R1", "9m_WT_R2")

df_select <- df %>% 
  filter(condition %in% compare_groups) %>% 
  mutate(genes = factor(genes, levels = ref_genes), condition = factor(condition, levels = compare_groups))

p <- ggplot(df_select, aes(x = genes, y = num_per_grid, fill = condition)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.75,
           colour = "#e9ecef",
           linewidth = 0.25) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Genes", y = "Count per Spot") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave(here::here("gene_capture_efficiency.jpeg"), p, width = 12, height = 5, dpi = dpi)

# Batch effect

df <- read.csv(here::here("count_per_gene_MERSCOPE_9m_WT.csv"))

cor_test_result <- cor.test(df$MERSCOPE_9m_WT_R1, df$MERSCOPE_9m_WT_R2, method = "pearson")
print(cor_test_result)

p <- ggplot(df, aes(x = MERSCOPE_9m_WT_R1, y = MERSCOPE_9m_WT_R2, label = Gene)) +
  geom_point(shape = 21, fill = "#a0ccec", color = "lightgrey", size = 1, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.75, linetype = "dashed", color = "darkblue") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = " ", x = " ", y = " ") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

ggsave(here::here("count_per_gene_MERSCOPE_9m_WT.jpeg"), p, width = 5, height = 5, dpi = dpi)


df <- read.csv(here::here("count_per_gene_MERSCOPE_9m_AD.csv"))

cor_test_result <- cor.test(df$MERSCOPE_9m_AD_R1, df$MERSCOPE_9m_AD_R2, method = "pearson")
print(cor_test_result)

p <- ggplot(df, aes(x = MERSCOPE_9m_AD_R1, y = MERSCOPE_9m_AD_R2, label = Gene)) +
  geom_point(shape = 21, fill = "#a0ccec", color = "lightgrey", size = 1, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.75, linetype = "dashed", color = "darkblue") +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = " ", x = " ", y = " ") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

ggsave(here::here("count_per_gene_MERSCOPE_9m_AD.jpeg"), p, width = 5, height = 5, dpi = dpi)



