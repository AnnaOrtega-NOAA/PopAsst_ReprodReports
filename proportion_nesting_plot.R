#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   ESA 2020 (NMFS and USFWS, 2020) 
#   Leatherback nesting over time, proportion each beach represents
#   Anna Ortega 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Set your working directory
setwd("C:/Users/anna.ortega/Desktop/PopAsst_RR")

library(tidyverse)
library(scales)

# 1. Load and Reshape
df <- read_csv("Dc_ESA_nestdat.csv")

df_long <- df %>%
  pivot_longer(
    cols = -Site, 
    names_to = "Year", 
    values_to = "Nests_raw",
    values_transform = list(Nests_raw = as.character)
  ) %>%
  mutate(
    Year = as.numeric(Year),
    Status = case_when(
      toupper(Nests_raw) == "X" ~ "Monitored (No Count)",
      !is.na(as.numeric(Nests_raw)) ~ "Count Recorded",
      TRUE ~ "No Data"
    ),
    Nests = as.numeric(Nests_raw)
  ) %>%
  filter(!is.na(Year), !is.na(Site))

# 2. Quantify Proportions
df_proportions <- df_long %>%
  group_by(Year) %>%
  mutate(Yearly_Total = sum(Nests, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Proportion = ifelse(Status == "Count Recorded" & Yearly_Total > 0, 
                        Nests / Yearly_Total, NA)
  )

# 3. Final Prep & REORDERING
df_final <- df_proportions %>%
  filter(Year == 1981 | Year >= 2002) %>%
  group_by(Site) %>%
  # Keep only sites that appear in our selected years
  filter(any(!is.na(Nests) | Status == "Monitored (No Count)")) %>%
  ungroup() %>%
  mutate(Site = factor(Site)) %>% # Reset factor levels
  group_by(Site) %>%
  # Calculate mean specifically for the years we are plotting
  mutate(avg_prop = mean(Proportion, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(avg_prop = replace_na(avg_prop, 0)) %>%
  # Rank: highest average proportion at the top
  mutate(Site = fct_reorder(Site, avg_prop)) %>% 
  mutate(
    Era = factor(ifelse(Year == 1981, "Past", "Modern"), levels = c("Past", "Modern")),
    Year_fact = factor(Year),
    label_text = case_when(
      Status == "Monitored (No Count)" ~ "X",
      !is.na(Proportion) & Proportion > 0 ~ percent(Proportion, accuracy = 1),
      TRUE ~ ""
    )
  )

# 4. Create the Visualization
p <- ggplot(df_final, aes(x = Year_fact, y = Site)) +
  # Background tiles
  geom_tile(fill = "white", color = "grey95", size = 0.1) +
  # Proportion tiles
  geom_tile(aes(fill = Proportion), color = "grey90", size = 0.2) +
  # Effort 'X' tiles
  geom_tile(data = filter(df_final, Status == "Monitored (No Count)"), 
            fill = "#FFD700", alpha = 0.8, color = "grey90") +
  # CELL LABELS (Increased size for readability)
  geom_text(aes(label = label_text), size = 3.2, fontface = "bold", 
            color = "black", family = "sans") +
  # Color Scale
  scale_fill_gradient(
    low = "white", 
    high = "dodgerblue4", 
    labels = percent_format(), 
    na.value = "transparent"
  ) +
  facet_grid(. ~ Era, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 12) + # Higher base size makes everything larger
  theme(
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    
    # Grid and Layout
    panel.grid = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    strip.text = element_blank(),
    axis.title = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
    
    # Legend at the bottom saves vertical space
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    plot.title = element_text(face = "bold", size = 14),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(
    fill = "Annual Share (%)",
    caption = "Yellow cells (X) indicate monitoring without counts. Source: ESA 2020 (NMFS & USFWS)"
  )

# 5. Display and Export
print(p)

# Save with large dimensions to ensure the text isn't cramped
ggsave("Leatherback_Nesting_Heatmap.png", p, 
       width = 16, height = 10, dpi = 300, bg = "white")