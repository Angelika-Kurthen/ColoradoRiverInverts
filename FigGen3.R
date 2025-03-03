###########
## Figure generation for 3rd manuscript
##########
library(ggpubr)
library(patchwork)
#install.packages("svglite")
library(svglite)
library(igraph)
library(ggraph)
library(dplyr)
library(broom) 
library(readr)

#source("Multispp.R")
source("MultisppHydropeak.R")

# make sure in the correct format
hydropeak.abunds$`rep(hydropeak[hyd], times = 5)` <- as.numeric(hydropeak.abunds$`rep(hydropeak[hyd], times = 5)`)
hydropeak.abunds$avg.abund<- as.numeric(hydropeak.abunds$avg.abund)
hydropeak.abunds$taxa <- as.factor(hydropeak.abunds$taxa)

#plot
hyd_abund_multi <- ggplot(data = hydropeak.abunds, aes(`rep(hydropeak[hyd], times = 5)`, avg.abund, group = taxa, color = taxa))+ 
  geom_line(linewidth = 1, alpha = 0.8)+
  # geom_ribbon(data = Hydro_Abund, aes(ymin = scale(means - sd),
  #                 ymax = scale(means + sd), fill = V4),
  #                        alpha = .15,
  #                        show.legend = T) +
  scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
  geom_vline(aes(xintercept = 0.01), linetype = "dotted", 
             linewidth=1)+
  geom_vline(aes(xintercept = 0.17 ), linetype="dashed", 
             linewidth=1)+
  geom_vline(aes(xintercept = 0.55 ), linetype = "dotdash", 
             linewidth=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme_bw()+
  xlab("Hydropeaking Intensity")+
  ylab("Relative Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

# make sure in correct format
hydropeak.biomass$`rep(hydropeak[hyd], times = 5)` <- as.numeric(hydropeak.biomass$`rep(hydropeak[hyd]`)
hydropeak.biomass$avg.biomass <- as.numeric(hydropeak.biomass$avg.biomass)
hydropeak.biomass$taxa <- as.factor(hydropeak.biomass$taxa) 
#plot
hyd_ts_multi <- ggplot(data = hydropeak.biomass, aes(`rep(hydropeak[hyd], times = 5)`,avg.biomass , color = taxa)) + 
  # geom_ribbon(aes(ymin = sizemeans - sizesd,
  #                   ymax = sizemeans + sizesd),
  #               colour = V3,
  #               alpha = .15,
  #               show.legend = T) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  xlab("Hydropeaking Intensity")+
  ylab("Relative Biomass (mg)")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

#make sure in correct format
# can also make this emergent by only looking at HYOS, CHIR, and BAET
hydropeak.avg.annual$`rep(hydropeak[hyd], times = 5)` <- as.numeric(hydropeak.avg.annual$`rep(hydropeak[hyd], times = 5)`)
hydropeak.avg.annual$taxa <- as.factor(hydropeak.avg.annual$taxa)
hydropeak.avg.annual$`means.s3.biomass$mean.S3.biomass` <- as.numeric(hydropeak.avg.annual$`means.s3.biomass$mean.S3.biomass`)
hyd_yr_multi <- ggplot(data = hydropeak.avg.annual, aes(`rep(hydropeak[hyd], times = 5)`, log(`means.s3.biomass$mean.S3.biomass`), group = taxa, color = taxa)) + 
  # geom_ribbon(aes(ymin = sizemeans - sizesd,
  #                   ymax = sizemeans + sizesd),
  #               colour = V3,
  #               alpha = .15,
  #               show.legend = T) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  xlab("Hydropeaking Intensity")+
  ylab("Log Annual S3 Biomass (mg)")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggarrange(hyd_abund_multi, hyd_ts_multi , hyd_yr_multi,
          labels = c("a", "b", "c"),
          ncol = 3, nrow = 1, common.legend = T)

# 8 different scenarios (proportions)

# read in the csv's
Multispp_temp_biomass <- read_csv("Multispp_temp_biomass.csv")
Multispp_temp_abund <- read_csv("Multispp_temp_abund(1).csv")
Multispp_temp_abund_spike <- read_csv("Multispp_temp_abund_spike.csv")
Multispp_temp_biomass_spike <- read_csv("Multispp_temp_biomass_spike.csv")
Multispp_temp_hyd_biomass <- read_csv("Multispp_temp_hyd_biomass.csv")
Multispp_temp_hyd_biomass_spike <- read_csv("Multispp_temp_hyd_biomass_spike.csv")
Multispp_temp_hyd_HFE_biomass_spike <- read_csv("Multispp_temp_hyd_HFE_biomass_spike.csv")
Multispp_temp_hyd_HFE_biomass <- read_csv("Multispp_temp_hyd_HFE_biomass.csv")
Multispp_temp_hyd_HFE_abund_spike <- read_csv("Multispp_temp_hyd_HFE_abund_spike.csv")
Multispp_temp_HFE_abund_spike <- read_csv("Multispp_temp_HFE_abund_spike.csv")
Multispp_temp_HFE_abund <- read_csv("Multispp_temp_HFE_abund.csv")
Multispp_temp_HFE_biomass <- read_csv("Multispp_temp_HFE_biomass.csv")
Multispp_temp_hyd_abund <- read_csv("Multispp_temp_hyd_abund.csv")
Multispp_temp_HFE_biomass_spike <- read_csv("Multispp_temp_HFE_biomass_spike.csv")
Multispp_temp_hyd_abund_spike <- read_csv("Multispp_temp_hyd_abund_spike.csv")
Multispp_temp_hyd_HFE_abund <- read_csv("Multispp_temp_hyd_HFE_abund.csv")

multi_abund_combo <- bind_rows(
  Multispp_temp_abund %>% mutate(source = "Temperature"),
  Multispp_temp_abund_spike  %>% mutate(source = "Temperature & spike"),
  Multispp_temp_hyd_abund %>% mutate(source = "Temperature & hydropeaking"),
  Multispp_temp_hyd_abund_spike %>% mutate(source = "Temperature & spike & hydropeaking"),
  #Multispp_temp_HFE_abund %>% mutate(source = "Temperature & HFE"),
  #Multispp_temp_HFE_abund_spike %>% mutate(source = "Temperature & spike & HFE"),
  Multispp_temp_hyd_HFE_abund  %>% mutate(source = "Temperature & hydropeaking & HFE"),
  Multispp_temp_hyd_HFE_abund_spike %>% mutate(source = "Temperature & spike & hydropeaking & HFE")
)


taxa.colors =c("#66CCEE", "#AA3377", "#228833", "#CCBB44", "#4477AA")

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  "Temperature +\nHydropeaking", 
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE", 
  "Temperature +\nSpike + HFE", 
  "Temperature +\nHydropeaking + HFE", 
  "Temperature +\nSpike + Hydropeaking + HFE"
)

multi_abund_combo$source <- ordered(multi_abund_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
multi_abund_combo$taxa <- ordered(multi_abund_combo$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))

system.time(
multi_abund_summary <- multi_abund_combo %>%
  group_by(temperature, taxa, source) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE),
            sd_abund = sd(abundance, na.rm = TRUE), .groups = "drop")
)
head(multi_abund_summary)
multi_abund <- ggplot(data = multi_abund_summary, 
                      aes(x = as.factor(temperature), y = mean_abund, color = taxa)) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(.~source, nrow = 4) +
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Abundance") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))

multi_biomass_combo <- bind_rows(
  Multispp_temp_biomass %>% mutate(source = "Temperature"),
  Multispp_temp_biomass_spike  %>% mutate(source = "Temperature & spike"),
  Multispp_temp_hyd_biomass %>% mutate(source = "Temperature & hydropeaking"),
  Multispp_temp_hyd_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking"),
  #Multispp_temp_HFE_biomass %>% mutate(source = "Temperature & HFE"),
  #Multispp_temp_HFE_biomass_spike %>% mutate(source = "Temperature & spike & HFE"),
  Multispp_temp_hyd_HFE_biomass  %>% mutate(source = "Temperature & hydropeaking & HFE"),
  Multispp_temp_hyd_HFE_biomass_spike %>% mutate(source = "Temperature & spike & hydropeaking & HFE")
)


label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  "Temperature +\nHydropeaking", 
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE", 
  "Temperature +\nSpike + HFE", 
  "Temperature +\nHydropeaking + HFE", 
  "Temperature +\nSpike + Hydropeaking + HFE"
)

multi_biomass_combo$source <- ordered(multi_biomass_combo$source, levels = c(
  "Temperature", 
  "Temperature & spike", 
  "Temperature & hydropeaking", 
  "Temperature & spike & hydropeaking",
  "Temperature & HFE", 
  "Temperature & spike & HFE", 
  "Temperature & hydropeaking & HFE", 
  "Temperature & spike & hydropeaking & HFE"))  
multi_biomass_combo$taxa <- ordered(multi_biomass_combo$taxa, levels = c("BAET", "HYOS", "CHIR", "GAMM", "NZMS"))

system.time(
  multi_biomass_summary <- multi_biomass_combo %>%
    group_by(temperature, taxa, source) %>%
    summarise(mean_biomass = mean(biomass, na.rm = TRUE),
              sd_biomass = sd(biomass, na.rm = TRUE), .groups = "drop")
)
head(multi_biomass_summary)
multi_biomass <- ggplot(data = multi_biomass_summary, 
                      aes(x = as.factor(temperature), y = mean_biomass, color = taxa)) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), 
                linewidth = 0.75, width = 0.2, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_discrete(labels=c("1" = "Baseline", "1.1" = "+10%",
                            "1.2" = "+20%", "1.5" = "+50%")) +
  facet_wrap(.~source, nrow= 4) +
  theme_bw() +
  scale_color_manual(name = " ", labels=c(expression(paste(italic("Baetidae")," spp.")), expression(paste(italic("Hydropsyche"), " spp.")),expression(paste(italic("Chironomidae"), " spp.")), expression(italic("G. lacustris")),  expression(italic("P. antipodarum"))), values=taxa.colors) +
  xlab("Temperature") +
  labs(y= "Relative Biomass") +
  theme(text = element_text(size = 13), 
        axis.text.x = element_text(hjust = 1, angle=45, size = 12), 
        axis.text.y = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent"))








Multispp_sensitivity_vitalrates <- read_csv("Multispp_sensitivity_vitalrates.csv")

sensgrid1 <- ggplot(Multispp_sensitivity_vitalrates, aes(x = SensitivityIncrement, y = Abundance, color = StageGroup)) +
  geom_point(size = 1) +  # Alternatively, you can add geom_line() for smoother trends.
  facet_grid(StageGroup ~ Parameter, scales = "free_y") +  # Facet by both StageGroup & Parameter
  theme_bw() +
  labs(
    x = "Sensitivity Increment",
    y = "Abundance",
    title = "Sensitivity Analysis by StageGroup and Parameter"
  ) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),  # Adjust facet label size
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"  # Hide legend (since color is redundant with facet_grid)
  )

ggsave("sensgrid1.png", sensgrid1, device = "png", width = 8, height = 11.5, units= c("in"), dpi = "retina")
# Fit a linear model for each StageGroup × Parameter combination
interaction_results <- Multispp_sensitivity_vitalrates %>%
  group_by(StageGroup, Parameter) %>%
  do({
    model <- lm(Abundance ~ SensitivityIncrement, data = .)  # Fit linear model
    model_summary <- summary(model)
    
    # Extract key values
    tibble(
      Slope = coef(model)[2],  # Extract slope
      P_value = coef(summary(model))["SensitivityIncrement", "Pr(>|t|)"],  # Get p-value for slope
      R2 = model_summary$r.squared  # Get R²
    )
  }) %>%
  ungroup() %>%
  # Filter results based on criteria
  filter(R2 > 0.3, P_value < 0.01)

# Print the filtered results
print(interaction_results)

library(igraph)
library(ggraph)
library(dplyr)
library(igraph)
library(ggraph)
library(dplyr)
library(stringr)


# Rename parameters: Change "G1" to "S1", "G2" to "S2"
interaction_results <- interaction_results %>%
  mutate(Parameter = case_when(
    grepl("G1_HYOS", Parameter) ~ "S1_HYOS",
    grepl("G2_HYOS", Parameter) ~ "S2_HYOS",
    grepl("G1_BAET", Parameter) ~ "S1_BAET",
    grepl("G2_BAET", Parameter) ~ "S2_BAET",
    grepl("G1_NZMS", Parameter) ~ "S1_NZMS",
    grepl("G2_NZMS", Parameter) ~ "S2_NZMS",
    grepl("G1_GAMM", Parameter) ~ "S1_GAMM",
    grepl("G2_GAMM", Parameter) ~ "S2_GAMM",
    grepl("G1_CHIR", Parameter) ~ "S1_CHIR",
    grepl("G2_CHIR", Parameter) ~ "S2_CHIR",
    TRUE ~ Parameter
  ))

# Extract the taxon type from parameter names
interaction_results <- interaction_results %>%
  mutate(Taxon = case_when(
    grepl("HYOS", Parameter) ~ "HYOS",
    grepl("GAMM", Parameter) ~ "GAMM",
    grepl("NZMS", Parameter) ~ "NZMS",
    grepl("CHIR", Parameter) ~ "CHIR",
    grepl("BAET", Parameter) ~ "BAET",
    TRUE ~ "Other"
  ))

# Create an edge list
edges <- interaction_results %>%
  select(Parameter, StageGroup, Slope, Taxon)

colnames(edges) <- c("from", "to", "weight", "taxon")  # Rename for igraph format

# Create a list of unique nodes (StageGroups & Parameters)
nodes <- data.frame(name = unique(c(edges$from, edges$to)))

# Extract taxon type for coloring
nodes <- nodes %>%
  mutate(
    Taxon = str_extract(name, "(HYOS|GAMM|NZMS|CHIR|BAET)"),
    Stage = str_extract(name, "S\\d+"),
    type = ifelse(is.na(Stage), "Parameter", "StageGroup"),
    label = ifelse(type == "StageGroup", Stage, name)  # Keep only "S1", "S2", "S3" for StageGroups
  )

# Ensure StageGroups are ordered correctly (S1 → S2 → S3)
nodes <- nodes %>%
  arrange(Taxon, Stage)

# Count total StageGroups
n_stagegroups <- sum(nodes$type == "StageGroup")

# Assign evenly spaced circular coordinates for StageGroups
nodes <- nodes %>%
  mutate(
    x = ifelse(type == "StageGroup", cos(seq(0, 2 * pi, length.out = n_stagegroups + 1))[1:n_stagegroups], NA),
    y = ifelse(type == "StageGroup", sin(seq(0, 2 * pi, length.out = n_stagegroups + 1))[1:n_stagegroups], NA)
  )

# Create igraph object
network_graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
# Define color palette for taxon-based nodes (matching arrows)
taxon_colors <- c(
  "HYOS" = "#AA3377",
  "GAMM" = "#CCBB44",
  "NZMS" = "#4477AA",
  "CHIR" = "#228833",
  "BAET" = "#66CCEE"
)


# Ensure `Taxon` exists in `network_layout`
network_layout <- create_layout(network_graph, layout = "nicely") %>%
  left_join(nodes, by = "name") %>%
  mutate(
    x = ifelse(!is.na(x.x), x.x, x.y),  # Use manually assigned x if available
    y = ifelse(!is.na(y.x), y.x, y.y)   # Use manually assigned y if available
  ) %>%
  select(name, x, y, Taxon, type, label, everything())  # Ensure `Taxon` is included
ggraph(network_layout) +
  geom_edge_arc(aes(edge_width = weight, color = taxon),
                arrow = arrow(length = unit(4, "mm")),
                end_cap = circle(4, "mm"),
                start_cap = circle(4, "mm"),
                alpha = 0.8,  # 80% transparency
                curvature = 0.2, 
                show.legend = F) +  #Curved edges
  geom_node_point(size = 8, aes(color = Taxon.x)) +  # Fix: Use `Taxon` from network_layout
  geom_node_text(aes(label = label.x), repel = F, size = 4, fontface = "bold") +  # Use short labels (S1, S2, S3)
  scale_edge_width(range = c(0.5, 3)) +  # Edge thickness based on R2
  scale_edge_color_manual(values = taxon_colors) +  # Taxon-based colors for edges
  scale_color_manual(labels =c("Baetidae spp.", "Chironomidae spp.", "G. lacustris", "Hydropsyche spp.", "P. antipodarum"), values = taxon_colors)+
  theme_void()+
  labs(title = " ", color = "Taxon", edge_width = "R² Strength")
