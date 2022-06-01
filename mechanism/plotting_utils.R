library(ggplot2)


add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


# data in tidy format
# columns: value, dataset/method, splice type (e.g.)

# modified from: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             # newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

# from: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


splittedViolinPlot = function(data) {
  
  data$db = factor(data$db, levels = c("ProteasomeDB", "randomDB"))
  
  # print stats
  s = data %>%
    dplyr::group_by(type, db) %>%
    dplyr::summarise(n = dplyr::n(),
                     mean = mean(value),
                     median = median(value),
                     std = sd(value))
  print.data.frame(s)
  
  theme_set(theme_classic())
  
  # plotting
  sc = ggplot(data = data, aes(x = type, y = value, fill = db)) +
    geom_split_violin(size = .25) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = .2,
                 position = position_dodge(width = .2)) +
    xlab("product type") +
    scale_fill_manual("database",
                      values = c(plottingCols[["ProteasomeDB"]], plottingCols[["randomDB"]]),
                      drop = F) +
    theme(axis.text.x = element_text(angle = 90))
  
  # add labs and title downstream
  
  return(sc)
}




splittedViolinPlot.dotProd = function(data) {
  
  # print stats
  s = data %>%
    group_by(type.DB, method) %>%
    summarise(n = dplyr::n(),
              mean = mean(value),
              median = median(value),
              std = sd(value))
  print.data.frame(s)
  
  theme_set(theme_classic())
  
  # plotting
  sc = ggplot(data = data, aes(x = type.DB, y = value, fill = method)) +
    geom_split_violin(size = .25) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = .2,
                 position = position_dodge(width = .2)) +
    xlab("product type") +
    scale_fill_manual("dataset and method",
                      values = c(plottingCols[["inSPIRE"]],plottingCols[["invitroSPI"]]),
                      drop = F) +
    theme(axis.text.x = element_text(angle = 90))
  
  # add labs and title downstream
  
  
  return(sc)
}
