# CLEAN DAG USING PURE GGPLOT2
# Full control over every element - no dagitty rendering issues
# Includes Wind-Rain correlation structure

library(ggplot2)
library(dplyr)
library(posterior)
library(scales)  # For alpha() function

# Load the fitted SCM model
mod.SEM.sim <- readRDS("Output/Models/mod.SEM.sim.rds")

# Extract posterior means for coefficients
coef_summary <- summarise_draws(mod.SEM.sim, "mean")

# Helper function to get coefficient value
get_coef <- function(name) {
  val <- coef_summary %>% 
    filter(variable == name) %>% 
    pull(mean)
  if(length(val) == 0) return(NA)
  return(round(val, 2))
}

# ============================================================================
# DEFINE NODES - Manual positioning
# ============================================================================

nodes <- data.frame(
  name = c("Rain", "Wind", "Substrate", "Micro", "Year", "Goby_lag",
           "U_phys", "U_bio",
           "BreachDays", "Temp", "DO", "SAV", "SC", "SB",
           "Goby"),
  x = c(1, 1, 1, 1, 1, 1,           # Exogenous (left)
        8.5, 9.5,                       # Latents (top right)
        3.5, 3.5, 3.5, 7, 7, 7,      # Mediators
        5.5),                          # Outcome (bottom right)
  y = c(8, 6, 4, 2.5, 1.5, 0.5,       # Exogenous
        9, 8,                           # Latents (top right: U_phys at 9, U_bio at 8)
        8, 6, 4, 1.5, 4, 5,      # Mediators
        0.5),                           # Outcome (bottom right)
  type = c(rep("Exogenous", 6),
           rep("Latent", 2),
           rep("Mediator", 6),
           "Outcome"),
  stringsAsFactors = FALSE
)

# ============================================================================
# DEFINE EDGES
# ============================================================================

edges <- data.frame(
  from = c("Rain", "Wind", "Wind", "BreachDays", "BreachDays", 
           "Temp", "Temp", "Temp",
           "DO", "DO", "DO", "DO", 
           "SAV", "SAV", "SAV", 
           "Substrate", "Substrate",
           "SC", "SB", 
           "Micro", "Year", "Goby_lag",
           "U_phys", "U_phys", "U_phys", 
           "U_bio", "U_bio", "U_bio"),
  to = c("BreachDays", "Temp", "DO", "Temp", "Goby", 
         "DO", "SAV", "Goby",
         "SAV", "SC", "SB", "Goby", 
         "SC", "SB", "Goby", 
         "SC", "Goby",
         "Goby", "Goby", 
         "Goby", "Goby", "Goby",
         "BreachDays", "Temp", "DO", 
         "SAV", "SC", "SB"),
  stringsAsFactors = FALSE
)

# Add coordinates
edges <- edges %>%
  left_join(nodes %>% select(name, x, y), by = c("from" = "name")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(nodes %>% select(name, x, y), by = c("to" = "name")) %>%
  rename(x2 = x, y2 = y)

# Shorten arrows to not overlap with nodes
node_radius <- 0.8
edges <- edges %>%
  mutate(
    dx = x2 - x1,
    dy = y2 - y1,
    len = sqrt(dx^2 + dy^2),
    x1_adj = x1 + (dx/len) * node_radius,
    y1_adj = y1 + (dy/len) * node_radius,
    x2_adj = x2 - (dx/len) * node_radius,
    y2_adj = y2 - (dy/len) * node_radius
  )

# Add Wind-Rain correlation (bidirectional arrow)
wind_rain_corr <- data.frame(
  x1 = 1, y1 = 6,
  x2 = 1, y2 = 8,
  x1_adj = 1, y1_adj = 6.35,
  x2_adj = 1, y2_adj = 7.65
)

# ============================================================================
# COLOR SCHEME
# ============================================================================

node_colors <- c(
  "Exogenous" = "#5B9BD5",  # Blue
  "Mediator" = "#70AD47",    # Green
  "Outcome" = "#C5504B",     # Red
  "Latent" = "#7F7F7F"       # Gray
)

node_fills <- c(
  "Exogenous" = "#5B9BD566",   # Blue with 40% alpha (66 in hex = ~40%)
  "Mediator" = "#70AD4766",     # Green with 40% alpha
  "Outcome" = "#C5504B66",      # Red with 40% alpha
  "Latent" = "#FFFFFF"          # White (fully opaque)
)

# ============================================================================
# CREATE PLOT
# ============================================================================

# Debug: Check nodes data
cat("\n=== Nodes Data ===\n")
print(nodes)
cat("\n=== Node Fills ===\n")
print(node_fills)

p <- ggplot() +
  # Draw causal arrows (BLACK) with open triangle heads
  geom_segment(
    data = edges,
    aes(x = x1_adj, y = y1_adj, xend = x2_adj, yend = y2_adj),
    arrow = arrow(length = unit(8, "pt"), type = "open", angle = 35),
    linewidth = 1.2,
    color = "black"
  ) +
  # Draw Wind-Rain correlation (double-headed arrow)
  geom_segment(
    data = wind_rain_corr,
    aes(x = x1_adj, y = y1_adj, xend = x2_adj, yend = y2_adj),
    arrow = arrow(length = unit(6, "pt"), type = "open", angle = 35, ends = "both"),
    linewidth = 0.8,
    color = "black",
    linetype = "dashed"
  ) +
  # Draw node circles - with alpha in color hex
  geom_point(
    data = nodes,
    aes(x = x, y = y, fill = type),
    size = 18,
    shape = 21,
    stroke = 2,
    color = "black"
  ) +
  # Add node text labels
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = name),
    size = 3.5,
    fontface = "plain",
    color = "black"
  ) +
  # Colors
  scale_fill_manual(
    name = "",
    values = node_fills,
    breaks = c("Exogenous", "Mediator", "Outcome", "Latent")
  ) +
  # Control legend symbol size
  guides(fill = guide_legend(override.aes = list(size = 9))) +  # Half of 18
  # Clean theme
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 10, margin = margin(b = 15), color = "gray40"),
    legend.position = "bottom",
    legend.text = element_text(size = 5),  # Reduced from 10 to 5 (50% smaller)
    plot.margin = margin(20, 20, 20, 20),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_fixed(ratio = 1) +
  labs(
    title = "Structural Causal Model: Tidewater Goby",
    subtitle = "Dashed line shows Wind-Rain correlation | Open circles are latent confounders"
  )

print(p)

# Save
ggsave("Output/Plots/dag_ggplot_clean.png", 
       plot = p, 
       width = 12, 
       height = 9, 
       dpi = 400,
       bg = "white")

# ============================================================================
# Version WITH coefficient labels
# ============================================================================

# Add coefficients to edges
edges$coef <- NA
edges$coef[edges$from == "Rain" & edges$to == "BreachDays"] <- get_coef("beta_Rain")
edges$coef[edges$from == "Wind" & edges$to == "Temp"] <- get_coef("beta_Wind_vec[2]")
edges$coef[edges$from == "Wind" & edges$to == "DO"] <- get_coef("beta_Wind_vec[1]")
edges$coef[edges$from == "Temp" & edges$to == "Goby"] <- get_coef("beta_Temp_vec[1]")
edges$coef[edges$from == "DO" & edges$to == "Goby"] <- get_coef("beta_DO_vec[1]")
edges$coef[edges$from == "SAV" & edges$to == "Goby"] <- get_coef("beta_SAV_vec[1]")
edges$coef[edges$from == "BreachDays" & edges$to == "Goby"] <- get_coef("beta_BreachDays_vec[1]")

# Calculate label positions
edges <- edges %>%
  mutate(
    x_mid = (x1 + x2) / 2,
    y_mid = (y1 + y2) / 2,
    label = ifelse(!is.na(coef), sprintf("%.2f", coef), "")
  )

p_labeled <- p +
  geom_label(
    data = edges %>% filter(label != ""),
    aes(x = x_mid, y = y_mid, label = label),
    size = 2.5,
    fill = "white",
    color = "black",
    alpha = 0.9,
    label.size = 0.15,
    label.padding = unit(0.12, "lines")
  )

print(p_labeled)

ggsave("Output/Plots/dag_ggplot_with_coefs.png", 
       plot = p_labeled, 
       width = 12, 
       height = 9, 
       dpi = 400,
       bg = "white")

# ============================================================================
# Summary
# ============================================================================

cat("\n=== DAG Plots Created ===\n")
cat("1. dag_ggplot_clean.png (structure only)\n")
cat("2. dag_ggplot_with_coefs.png (with key coefficients)\n\n")

cat("=== Key Features ===\n")
cat("- BLACK arrows (linewidth = 1.2, color = black)\n")
cat("- Latent U's are OPEN circles (white fill, gray border)\n")
cat("- Wind-Rain correlation shown as dashed double-headed arrow\n")
cat("- Text labels with colored circle backgrounds\n")
cat("- Clean ggplot2 - no rendering issues\n\n")

cat("=== Complete! ===\n")