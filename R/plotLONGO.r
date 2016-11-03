#' main function
#'
#' analyzes data by creating windows and moving by a step size
#' @param simp.df dataframe to be plotted
plotLONGO <- function(simp.df) {
  x_vals <- (simp.df$kb_length)
  ymax <- (max(simp.df[, 2:(ncol(simp.df))]) * 1.2)
  yminim <- min(simp.df[, 2:(ncol(simp.df))])
  png(
    "LONGO_out.png",
    width = 6,
    height = 6,
    units = "in",
    res = 300
  )
  matplot(
    x = x_vals,
    y = simp.df[, 2],
    type = "l",
    col = 1,
    xlab = "gene length",
    ylab = "expression",
    ylim = c(yminim, ymax),
    main = "LONGO output"
  )
  for (i in 3:ncol(simp.df)) {
    par(new = TRUE)
    matplot(
      x = x_vals,
      y = simp.df[, i],
      type = "l",
      col = i - 1,
      xlab = "",
      ylab = "",
      ylim = c(yminim, ymax),
      axes = FALSE
    )
  }
  labels <- colnames(simp.df)
  legend(
    "topright",
    legend = c(labels[2:length(labels)]),
    col = 2:ncol(simp.df) - 1,
    lty = 1,
    cex = 0.8,
    ncol = 2
  )
  dev.off()
}
