# ==============================================================================
# Test Publication-Style Recovery Plots
# ==============================================================================
# Generates recovery plots matching Haines et al. (2018) ORL paper style.
#
# Usage: Rscript analysis/2_plotting/test_recovery_plots.R
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, tidyr, patchwork)

source("analysis/2_plotting/plotting_publication.R")

# ==============================================================================
# Load Pseudo Data
# ==============================================================================
output_dir <- "outputs/recovery"

df_pvl <- read.csv(file.path(output_dir, "pseudo_recovery_pvl_delta.csv"))
df_orl <- read.csv(file.path(output_dir, "pseudo_recovery_orl.csv"))
df_eef <- read.csv(file.path(output_dir, "pseudo_recovery_eef.csv"))

cat("Loaded pseudo data for plotting.\n")

# ==============================================================================
# Clear existing plots
# ==============================================================================
plot_dir <- "plots/recovery"
if (dir.exists(plot_dir)) {
    unlink(plot_dir, recursive = TRUE)
    cat("Cleared existing plot folder.\n")
}
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Convert to long format
# ==============================================================================
df_pvl_long <- recovery_wide_to_long(df_pvl)
df_orl_long <- recovery_wide_to_long(df_orl)
df_eef_long <- recovery_wide_to_long(df_eef)

# ==============================================================================
# Generate Publication-Style Plots
# ==============================================================================
cat("Generating PVL-Delta recovery plot (Haines style)...\n")
p_pvl <- plot_recovery_publication(df_pvl_long, "PVL-Delta")

cat("Generating ORL recovery plot (Haines style)...\n")
p_orl <- plot_recovery_publication(df_orl_long, "ORL")

cat("Generating EEF recovery plot (Haines style)...\n")
p_eef <- plot_recovery_publication(df_eef_long, "EEF")

# ==============================================================================
# Save Plots
# ==============================================================================
plot_dir <- "plots/recovery"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(plot_dir, "recovery_pvl_delta_publication.png"), p_pvl,
    width = 12, height = 6, dpi = 200
)
ggsave(file.path(plot_dir, "recovery_orl_publication.png"), p_orl,
    width = 12, height = 6, dpi = 200
)
ggsave(file.path(plot_dir, "recovery_eef_publication.png"), p_eef,
    width = 12, height = 6, dpi = 200
)

# Also save as PDF for publication
ggsave(file.path(plot_dir, "recovery_pvl_delta_publication.pdf"), p_pvl,
    width = 12, height = 6
)
ggsave(file.path(plot_dir, "recovery_orl_publication.pdf"), p_orl,
    width = 12, height = 6
)
ggsave(file.path(plot_dir, "recovery_eef_publication.pdf"), p_eef,
    width = 12, height = 6
)

cat("\n=== Publication-Style Recovery Plots Saved ===\n")
cat("Output directory:", plot_dir, "\n")
cat("PNG files:\n")
cat("  - recovery_pvl_delta_publication.png\n")
cat("  - recovery_orl_publication.png\n")
cat("  - recovery_eef_publication.png\n")
cat("PDF files:\n")
cat("  - recovery_pvl_delta_publication.pdf\n")
cat("  - recovery_orl_publication.pdf\n")
cat("  - recovery_eef_publication.pdf\n")

# ==============================================================================
# Create Combined Multi-Model Figure (Like Fig. 4 in Haines)
# ==============================================================================
cat("\nGenerating combined figure...\n")

# Stack all models vertically
combined_all <- p_pvl / p_orl / p_eef +
    plot_annotation(
        title = "Parameter Recovery Results",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )

ggsave(file.path(plot_dir, "recovery_all_models_publication.png"), combined_all,
    width = 12, height = 14, dpi = 200
)
ggsave(file.path(plot_dir, "recovery_all_models_publication.pdf"), combined_all,
    width = 12, height = 14
)

cat("  - recovery_all_models_publication.png\n")
cat("  - recovery_all_models_publication.pdf\n")
