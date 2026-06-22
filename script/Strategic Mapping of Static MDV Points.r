# ===============================
# 1. LOAD LIBRARIES
# ===============================

library(sf)
library(terra)
library(dplyr)
library(exactextractr)
library(readr)
library(ggplot2)
library(ggspatial)
library(tidyr)
library(viridis)
library(scales)

# ===============================
# 2. IMPORT DATA
# ===============================

schools <- st_read("~/Documents/RABIES_ANALYSIS/inputs/cleaned_schools.gpkg")

pop_rast <- rast("~/Documents/RABIES_ANALYSIS/inputs/dog_population.tif")

urban_rast <- rast("~/Documents/RABIES_ANALYSIS/inputs/DEGREE OF URBANISATION 2026 RECLASSIFIED.tif")

country <- st_read("~/Documents/RABIES_ANALYSIS/inputs/Kenya.gpkg")

counties <- st_read("~/Documents/RABIES_ANALYSIS/inputs/counties.gpkg")

subcounties <- st_read("~/Documents/RABIES_ANALYSIS/inputs/Kenya_subcounties.gpkg")

output_dir <- "~/Documents/RABIES_ANALYSIS/outputs/paper 3 outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ===============================
# 3. SETTINGS
# ===============================

target_crs <- "EPSG:32737"
map_crs    <- "EPSG:4326"

county_name_col    <- "county"
subcounty_name_col <- "subcounty"

buffer_scenarios <- data.frame(
  scenario = c("1km", "3km"),
  distance_m = c(1000, 3000)
)

# ===============================
# 4. REPROJECT DATA
# ===============================

schools     <- st_transform(schools, target_crs)
country     <- st_transform(country, target_crs)
counties    <- st_transform(counties, target_crs)
subcounties <- st_transform(subcounties, target_crs)

pop_rast <- project(pop_rast, target_crs, method = "near")

urban_rast <- project(urban_rast, pop_rast, method = "near")

# ===============================
# 5. BASIC CHECKS
# ===============================

if (!(county_name_col %in% names(counties))) {
  stop("County name column not found: ", county_name_col)
}

if (!(subcounty_name_col %in% names(subcounties))) {
  stop("Subcounty name column not found: ", subcounty_name_col)
}

total_pop <- terra::global(pop_rast, "sum", na.rm = TRUE)[1, 1]
cat("Total dog population:", total_pop, "\n")

# ===============================
# 6. HELPER FUNCTIONS
# ===============================

safe_exact_sum <- function(raster_obj, polygon_obj) {
  if (is.null(polygon_obj) || nrow(polygon_obj) == 0) return(0)
  
  val <- exactextractr::exact_extract(raster_obj, polygon_obj, "sum")
  val <- sum(as.numeric(val), na.rm = TRUE)
  
  ifelse(is.na(val), 0, val)
}

safe_exact_count <- function(raster_obj, polygon_obj, class_value) {
  if (is.null(polygon_obj) || nrow(polygon_obj) == 0) return(0)
  
  class_rast <- raster_obj == class_value
  
  val <- exactextractr::exact_extract(class_rast, polygon_obj, "sum")
  val <- sum(as.numeric(val), na.rm = TRUE)
  
  ifelse(is.na(val), 0, val)
}

save_map <- function(plot_obj, filename, width = 10, height = 8) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

make_label_points <- function(sf_obj) {
  suppressWarnings(st_point_on_surface(sf_obj))
}

# ===============================
# 7. CREATE RURAL AND URBAN DOG POPULATION RASTERS
# ===============================

rural_pop_rast <- ifel(urban_rast == 1, pop_rast, NA)
urban_pop_rast <- ifel(urban_rast == 2, pop_rast, NA)

# ===============================
# 8. MAP DATA
# ===============================

schools_map     <- st_transform(schools, map_crs)
counties_map    <- st_transform(counties, map_crs)
subcounties_map <- st_transform(subcounties, map_crs)

county_labels_map <- make_label_points(counties_map)

# ===============================
# 9. MAP 1: SCHOOLS AND COUNTIES
# ===============================

p_schools_counties <- ggplot() +
  geom_sf(data = counties_map, fill = NA, color = "grey40", linewidth = 0.3) +
  geom_sf(data = schools_map, aes(color = "Schools"), size = 0.6, alpha = 0.8) +
  geom_sf_text(
    data = county_labels_map,
    aes(label = .data[[county_name_col]]),
    size = 1.8,
    color = "black",
    check_overlap = TRUE
  ) +
  scale_color_manual(values = c("Schools" = "blue"), name = "Legend") +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +
  labs(title = "Kenya Schools") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  )

save_map(p_schools_counties, "map_schools_and_counties.png")

# ===============================
# 10. ADMIN ANALYSIS FUNCTION
# ===============================

run_admin_analysis <- function(admin_sf, admin_name_col, admin_level, buffer_dissolved,
                               scenario_name, buffer_dist, pop_buffer) {
  
  admin_results <- list()
  
  for (i in seq_len(nrow(admin_sf))) {
    
    admin_area <- admin_sf[i, ]
    admin_name <- admin_area[[admin_name_col]][1]
    
    admin_total_pop <- safe_exact_sum(pop_rast, admin_area)
    
    buffer_clip <- suppressWarnings(st_intersection(buffer_dissolved, admin_area))
    buffer_clip <- st_make_valid(buffer_clip)
    
    admin_covered_pop <- safe_exact_sum(pop_rast, buffer_clip)
    admin_uncovered_pop <- admin_total_pop - admin_covered_pop
    
    uncovered_area <- suppressWarnings(st_difference(admin_area, st_union(buffer_dissolved)))
    uncovered_area <- st_make_valid(uncovered_area)
    
    uncovered_rural_pop <- safe_exact_sum(rural_pop_rast, uncovered_area)
    uncovered_urban_pop <- safe_exact_sum(urban_pop_rast, uncovered_area)
    
    uncovered_rural_pixel_count <- safe_exact_count(urban_rast, uncovered_area, 1)
    uncovered_urban_pixel_count <- safe_exact_count(urban_rast, uncovered_area, 2)
    
    coverage_pct <- ifelse(
      admin_total_pop > 0,
      (admin_covered_pop / admin_total_pop) * 100,
      NA
    )
    
    admin_results[[i]] <- data.frame(
      scenario = scenario_name,
      buffer_distance_m = buffer_dist,
      admin_level = admin_level,
      admin_name = admin_name,
      total_dog_population = admin_total_pop,
      covered_dog_population = admin_covered_pop,
      uncovered_dog_population = admin_uncovered_pop,
      coverage_pct = coverage_pct,
      pct_of_country_dog_population = (admin_total_pop / total_pop) * 100,
      pct_of_total_buffer_population = ifelse(
        pop_buffer > 0,
        (admin_covered_pop / pop_buffer) * 100,
        NA
      ),
      
      uncovered_rural_dog_population = uncovered_rural_pop,
      uncovered_urban_dog_population = uncovered_urban_pop,
      
      pct_uncovered_population_rural = ifelse(
        admin_uncovered_pop > 0,
        (uncovered_rural_pop / admin_uncovered_pop) * 100,
        NA
      ),
      
      pct_uncovered_population_urban = ifelse(
        admin_uncovered_pop > 0,
        (uncovered_urban_pop / admin_uncovered_pop) * 100,
        NA
      ),
      
      uncovered_rural_pixel_count = uncovered_rural_pixel_count,
      uncovered_urban_pixel_count = uncovered_urban_pixel_count
    )
  }
  
  bind_rows(admin_results) %>%
    mutate(
      coverage_class = case_when(
        coverage_pct >= 70 ~ "≥70% coverage",
        coverage_pct >= 50 & coverage_pct < 70 ~ "50–69% coverage",
        coverage_pct < 50 ~ "<50% coverage",
        TRUE ~ NA_character_
      )
    )
}

# ===============================
# 11. RUN BUFFER SCENARIOS
# ===============================

all_county_results <- list()
all_subcounty_results <- list()
all_buffers <- list()

for (s in seq_len(nrow(buffer_scenarios))) {
  
  scenario_name <- buffer_scenarios$scenario[s]
  buffer_dist   <- buffer_scenarios$distance_m[s]
  
  cat("\nRunning scenario:", scenario_name, "\n")
  
  buffers <- st_buffer(schools, dist = buffer_dist)
  
  buffer_dissolved <- st_union(buffers) %>%
    st_as_sf() %>%
    st_set_crs(target_crs) %>%
    st_make_valid()
  
  all_buffers[[scenario_name]] <- buffer_dissolved
  
  pop_buffer <- safe_exact_sum(pop_rast, buffer_dissolved)
  pop_gap <- total_pop - pop_buffer
  
  cat("Population inside", scenario_name, "buffer:", pop_buffer, "\n")
  cat("Population outside", scenario_name, "buffer:", pop_gap, "\n")
  
  # ===============================
  # COUNTY ANALYSIS
  # ===============================
  
  county_df <- run_admin_analysis(
    admin_sf = counties,
    admin_name_col = county_name_col,
    admin_level = "county",
    buffer_dissolved = buffer_dissolved,
    scenario_name = scenario_name,
    buffer_dist = buffer_dist,
    pop_buffer = pop_buffer
  ) %>%
    rename(county = admin_name)
  
  write_csv(
    county_df,
    file.path(output_dir, paste0(scenario_name, "_county_population_coverage.csv"))
  )
  
  all_county_results[[scenario_name]] <- county_df
  
  # ===============================
  # SUBCOUNTY ANALYSIS
  # ===============================
  
  subcounty_df <- run_admin_analysis(
    admin_sf = subcounties,
    admin_name_col = subcounty_name_col,
    admin_level = "subcounty",
    buffer_dissolved = buffer_dissolved,
    scenario_name = scenario_name,
    buffer_dist = buffer_dist,
    pop_buffer = pop_buffer
  ) %>%
    rename(subcounty = admin_name)
  
  write_csv(
    subcounty_df,
    file.path(output_dir, paste0(scenario_name, "_subcounty_population_coverage.csv"))
  )
  
  all_subcounty_results[[scenario_name]] <- subcounty_df
  
  # ===============================
  # MAP 2: BUFFER + COUNTIES
  # ===============================
  
  buffer_map <- st_transform(buffer_dissolved, map_crs)
  
  p_pop_buffer <- ggplot() +
    geom_sf(
      data = buffer_map,
      aes(color = "Dissolved buffer"),
      fill = "blue",
      alpha = 0.5,
      linewidth = 0.2
    ) +
    geom_sf(
      data = counties_map,
      aes(color = "Counties"),
      fill = NA,
      linewidth = 0.3
    ) +
    geom_sf_text(
      data = county_labels_map,
      aes(label = .data[[county_name_col]]),
      size = 1.8,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_color_manual(
      values = c(
        "Dissolved buffer" = "blue",
        "Counties" = "black"
      ),
      name = "Legend"
    ) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = north_arrow_fancy_orienteering
    ) +
    labs(title = paste0("Counties and ", scenario_name, " School Buffer")) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  save_map(
    p_pop_buffer,
    paste0("buffer_and_counties_", scenario_name, ".png")
  )
  
  # ===============================
  # MAP 3: COUNTY COVERAGE CLASS
  # ===============================
  
  counties_scenario_map <- counties_map %>%
    left_join(
      county_df,
      by = setNames("county", county_name_col)
    )
  
  county_labels <- make_label_points(counties_scenario_map)
  
  p_coverage_class <- ggplot() +
    geom_sf(
      data = counties_scenario_map,
      aes(fill = coverage_class),
      color = "grey30",
      linewidth = 0.3
    ) +
    geom_sf_text(
      data = county_labels,
      aes(label = .data[[county_name_col]]),
      size = 1.8,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_fill_manual(
      values = c(
        "≥70% coverage" = "green",
        "50–69% coverage" = "orange",
        "<50% coverage" = "red"
      ),
      name = "Coverage class",
      drop = FALSE
    ) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = north_arrow_fancy_orienteering
    ) +
    labs(
      title = paste0("County-Level Dog Population Coverage: ", scenario_name, " Buffer")
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  save_map(
    p_coverage_class,
    paste0("county_coverage_class_", scenario_name, ".png")
  )
  
  # ===============================
  # MAP 4: % UNCOVERED DOG POPULATION BY COUNTY
  # ===============================
  
  counties_scenario_map <- counties_scenario_map %>%
    mutate(
      uncovered_pct = ifelse(
        total_dog_population > 0,
        (uncovered_dog_population / total_dog_population) * 100,
        NA
      )
    )
  
  county_labels <- make_label_points(counties_scenario_map)
  
  p_uncovered <- ggplot() +
    geom_sf(
      data = counties_scenario_map,
      aes(fill = uncovered_pct),
      color = "grey30",
      linewidth = 0.3
    ) +
    geom_sf_text(
      data = county_labels,
      aes(label = .data[[county_name_col]]),
      size = 1.8,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_fill_gradient(
      low = "green",
      high = "red",
      labels = function(x) paste0(round(x, 1), "%"),
      name = "% Uncovered\ndog population"
    ) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = north_arrow_fancy_orienteering
    ) +
    labs(
      title = paste0(
        "Percentage of Uncovered Dog Population by County: ",
        scenario_name,
        " Buffer"
      )
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  save_map(
    p_uncovered,
    paste0("uncovered_dog_population_pct_", scenario_name, ".png")
  )
  
  # ===============================
  # MAP 5: SUBCOUNTY COVERAGE CLASS
  # ===============================
  
  subcounties_scenario_map <- subcounties_map %>%
    left_join(
      subcounty_df,
      by = setNames("subcounty", subcounty_name_col)
    )
  
  subcounty_labels <- make_label_points(subcounties_scenario_map)
  
  p_subcounty_coverage <- ggplot() +
    geom_sf(
      data = subcounties_scenario_map,
      aes(fill = coverage_class),
      color = "grey30",
      linewidth = 0.15
    ) +
    geom_sf_text(
      data = subcounty_labels,
      aes(label = .data[[subcounty_name_col]]),
      size = 1.8,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_fill_manual(
      values = c(
        "≥70% coverage" = "green",
        "50–69% coverage" = "orange",
        "<50% coverage" = "red"
      ),
      name = "Coverage class",
      drop = FALSE
    ) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = north_arrow_fancy_orienteering
    ) +
    labs(
      title = paste0("Subcounty-Level Dog Population Coverage: ", scenario_name, " Buffer")
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  save_map(
    p_subcounty_coverage,
    paste0("subcounty_coverage_class_", scenario_name, ".png"),
    width = 12,
    height = 9
  )
  
  # ===============================
  # MAP 6: % UNCOVERED DOG POPULATION BY SUBCOUNTY
  # ===============================
  
  subcounties_scenario_map <- subcounties_scenario_map %>%
    mutate(
      uncovered_pct = ifelse(
        total_dog_population > 0,
        (uncovered_dog_population / total_dog_population) * 100,
        NA
      )
    )
  
  subcounty_labels <- make_label_points(subcounties_scenario_map)
  
  p_subcounty_uncovered <- ggplot() +
    geom_sf(
      data = subcounties_scenario_map,
      aes(fill = uncovered_pct),
      color = "grey30",
      linewidth = 0.15
    ) +
    geom_sf_text(
      data = subcounty_labels,
      aes(label = .data[[subcounty_name_col]]),
      size = 1.8,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_fill_gradient(
      low = "green",
      high = "red",
      labels = function(x) paste0(round(x, 1), "%"),
      name = "% Uncovered\ndog population"
    ) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = north_arrow_fancy_orienteering
    ) +
    labs(
      title = paste0(
        "Percentage of Uncovered Dog Population by Subcounty: ",
        scenario_name,
        " Buffer"
      )
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  save_map(
    p_subcounty_uncovered,
    paste0("subcounty_uncovered_dog_population_pct_", scenario_name, ".png"),
    width = 12,
    height = 9
  )
}

# ===============================
# 12. COMBINE AND EXPORT RESULTS
# ===============================

final_county_df <- bind_rows(all_county_results)

write_csv(
  final_county_df,
  file.path(output_dir, "combined_1km_3km_county_population_coverage.csv")
)

final_subcounty_df <- bind_rows(all_subcounty_results)

write_csv(
  final_subcounty_df,
  file.path(output_dir, "combined_1km_3km_subcounty_population_coverage.csv")
)

# ===============================
# 13. BAR CHART: COUNTY COVERAGE %, 1KM VS 3KM
# ===============================

bar_df <- final_county_df %>%
  mutate(
    county = factor(county, levels = unique(county[order(coverage_pct)]))
  )

bar_colors <- c(
  "1km" = "#1f78b4",
  "3km" = "#33a02c"
)

p_bar <- ggplot(
  bar_df,
  aes(x = county, y = coverage_pct, fill = scenario)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(
    yintercept = 70,
    linetype = "dashed",
    linewidth = 0.8,
    color = "red"
  ) +
  annotate(
    "text",
    x = 1,
    y = 72,
    label = "70% target",
    hjust = 0,
    color = "red",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = bar_colors,
    name = "Buffer scenario"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, max(100, max(bar_df$coverage_pct, na.rm = TRUE)))
  ) +
  labs(
    title = "County Dog Population Coverage by School Buffer Distance",
    x = "County",
    y = "Coverage (%)"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

save_map(
  p_bar,
  "bar_chart_county_coverage_1km_vs_3km.png",
  width = 16,
  height = 8
)

# ===============================
# 14. BAR CHART: SUBCOUNTY COVERAGE %, 1KM VS 3KM
# ===============================

subcounty_bar_df <- final_subcounty_df %>%
  mutate(
    subcounty = factor(subcounty, levels = unique(subcounty[order(coverage_pct)]))
  )

p_subcounty_bar <- ggplot(
  subcounty_bar_df,
  aes(x = subcounty, y = coverage_pct, fill = scenario)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(
    yintercept = 70,
    linetype = "dashed",
    linewidth = 0.8,
    color = "red"
  ) +
  scale_fill_manual(
    values = bar_colors,
    name = "Buffer scenario"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, max(100, max(subcounty_bar_df$coverage_pct, na.rm = TRUE)))
  ) +
  labs(
    title = "Subcounty Dog Population Coverage by School Buffer Distance",
    x = "Subcounty",
    y = "Coverage (%)"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
  )

save_map(
  p_subcounty_bar,
  "bar_chart_subcounty_coverage_1km_vs_3km.png",
  width = 24,
  height = 9
)

# ===============================
# 15. SAVE DISSOLVED BUFFERS
# ===============================

for (scenario_name in names(all_buffers)) {
  
  buffer_out <- st_transform(all_buffers[[scenario_name]], map_crs)
  
  st_write(
    buffer_out,
    file.path(output_dir, paste0("dissolved_school_buffer_", scenario_name, ".gpkg")),
    delete_dsn = TRUE,
    quiet = TRUE
  )
}

# ===============================
# 16. FINAL MESSAGE
# ===============================

cat("\n✅ Analysis complete.\n")
cat("Outputs saved in:\n", output_dir, "\n")
cat("\nMain outputs:\n")
cat(" - combined_1km_3km_county_population_coverage.csv\n")
cat(" - combined_1km_3km_subcounty_population_coverage.csv\n")
cat(" - 1km_county_population_coverage.csv\n")
cat(" - 3km_county_population_coverage.csv\n")
cat(" - 1km_subcounty_population_coverage.csv\n")
cat(" - 3km_subcounty_population_coverage.csv\n")
cat(" - map_schools_and_counties.png\n")
cat(" - buffer_and_counties_1km.png\n")
cat(" - buffer_and_counties_3km.png\n")
cat(" - county_coverage_class_1km.png\n")
cat(" - county_coverage_class_3km.png\n")
cat(" - uncovered_dog_population_pct_1km.png\n")
cat(" - uncovered_dog_population_pct_3km.png\n")
cat(" - subcounty_coverage_class_1km.png\n")
cat(" - subcounty_coverage_class_3km.png\n")
cat(" - subcounty_uncovered_dog_population_pct_1km.png\n")
cat(" - subcounty_uncovered_dog_population_pct_3km.png\n")
cat(" - bar_chart_county_coverage_1km_vs_3km.png\n")
cat(" - bar_chart_subcounty_coverage_1km_vs_3km.png\n")