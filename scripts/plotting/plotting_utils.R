# shared plotting functions for recovery scripts

library(ggplot2)

# scatter plot with identity line for comparing true vs recovered params
plot_recovery <- function(true_vals, infer_vals, param_name) {
    df <- data.frame(true = true_vals, inferred = infer_vals)
    r_val <- round(cor(true_vals, infer_vals, use = "complete.obs"), 3)
    rmse_val <- round(sqrt(mean((true_vals - infer_vals)^2, na.rm = TRUE)), 3)

    ggplot(df, aes(x = true, y = inferred)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.5) +
        labs(
            title = param_name,
            subtitle = paste0("r = ", r_val, ", RMSE = ", rmse_val),
            x = "True Value",
            y = "Recovered Value"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
}
