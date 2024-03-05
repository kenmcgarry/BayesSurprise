
spacka <- function (seqdata, indic = "cplx", window.size = 0.2, sliding = TRUE, 
          wstep = 1, with.missing = FALSE, endmiss.as.void = FALSE, 
          silent.indic = TRUE, ...) 
{
  if (!inherits(seqdata, "stslist")) 
    TraMineR:::msg.stop("data is NOT a sequence object, see seqdef function to create one")
  if (length(indic) > 1) {
    if (indic[1] != "Develop") {
      TraMineR:::msg.warn("Vector indic, only first value is used!")
      indic <- iindic <- indic[1]
    }
    else if (length(indic) > 2) {
      TraMineR:::msg.warn("Develop vector indic, only first two value are used!")
      indic <- indic[1:2]
      iiindic <- iindic[2]
    }
    else iindic <- indic[2]
  }
  else iindic <- indic
  group.list <- c("all", "basic", "diversity", "complexity", 
                  "binary", "ranked")
  if (iindic %in% group.list) 
    TraMineR:::msg.stop("Bad indic value, group name not supported!")
  wstep <- as.integer(wstep)
  if (wstep < 1) 
    TraMineR:::msg.stop("wstep must be a strictly positive integer!")
  lgth <- seqlength(seqdata, with.missing = TRUE)
  maxl <- max(lgth)
  if (window.size <= 0) 
    TraMineR:::msg.stop("bad window.size value!")
  if (window.size < 1) 
    step <- ceiling(window.size * maxl)
  else {
    step <- min(window.size, maxl)
  }
  re <- step
  rs <- 1
  slid <- as.integer(sliding)
  nr <- nrow(seqdata)
  nc <- ceiling((maxl - step + 1)/wstep)
  j <- 1
  if (sliding) {
    if (silent.indic) {
      windic <- function(k, seqdata, indic, with.missing, 
                         ...) {
        suppressMessages(seqindic(seqdata[, (k - step + 
                                               1):k], indic = indic, with.missing = with.missing, 
                                  ...))
      }
    }
    else {
      windic <- function(k, seqdata, indic, with.missing, 
                         ...) {
        seqindic(seqdata[, (k - step + 1):k], indic = indic, 
                 with.missing = with.missing, ...)
      }
    }
  }
  else {
    if (silent.indic) {
      windic <- function(k, seqdata, indic, with.missing, 
                         ...) {
        suppressMessages(seqindic(seqdata[, 1:k], indic = indic, 
                                  with.missing = with.missing, ...))
      }
    }
    else {
      windic <- function(k, seqdata, indic, with.missing, 
                         ...) {
        seqindic(seqdata[, 1:k], indic = indic, with.missing = with.missing, 
                 ...)
      }
    }
  }
  ind.dyn <- matrix(NA, nrow = nr, ncol = nc)
  re.range <- seq(from = re, to = maxl, by = wstep)
  if (length(re.range) < 2) 
    TraMineR:::msg.stop("There is only one window, wstep probably too large!")
  void <- attr(seqdata, "void")
  nr <- attr(seqdata, "nr")
  miss.code <- if (endmiss.as.void & !with.missing) 
    c(void, nr)
  else void
  j <- 0
  for (k in re.range) {
    j <- j + 1
    ind.dyn[, j] <- as.matrix(windic(k, seqdata = seqdata, 
                                     indic = indic, with.missing = with.missing, ...))
    wna <- seqdata[, k] %in% miss.code
    ind.dyn[wna, j] <- NA
  }
  colnames(ind.dyn) <- colnames(seqdata)[re.range]
  rownames(ind.dyn) <- rownames(seqdata)
  class(ind.dyn) <- c("dynin", class(ind.dyn))
  attr(ind.dyn, "weights") <- attr(seqdata, "weights")
  attr(ind.dyn, "xtstep") <- attr(seqdata, "xtstep")
  attr(ind.dyn, "tick.last") <- attr(seqdata, "tick.last")
  attr(ind.dyn, "window.size") <- step
  attr(ind.dyn, "wstep") <- wstep
  attr(ind.dyn, "sliding") <- sliding
  attr(ind.dyn, "indic") <- iindic
  return(ind.dyn)}


####################

bs_plot <- function (seqdata = NULL, diss = NULL, k = NULL, sortv = "mds", 
          weighted = TRUE, grp.meth = "prop", squared = FALSE, pow = NULL, 
          seqrfobject = NULL, border = FALSE, ylab = NULL, yaxis = TRUE, 
          which.plot = "both", quality = TRUE, box.color = NULL, box.fill = NULL, 
          box.alpha = NULL, outlier.jitter.height = 0, outlier.color = NULL, 
          outlier.fill = NULL, outlier.shape = 19, outlier.size = 1.5, 
          outlier.stroke = 0.5, outlier.alpha = NULL) 
{
  if (inherits(seqrfobject, "seqrf") & inherits(seqdata, "stslist")) {
    usethis::ui_info("you specified a {usethis::ui_code('seqrfobject')} & {usethis::ui_code('seqdata')};\n    the latter as well as the potentially specified parameters\n    {usethis::ui_code(c('k', 'sortv', 'weighted', 'grp.meth', 'squared', 'pow'))} will be ignored;\n    the plot will be rendered for the {usethis::ui_field('seqrfobject')}")
  }
  if (!is.null(seqdata) & !inherits(seqdata, "stslist") & !inherits(seqdata, 
                                                                    "seqrf") & inherits(seqrfobject, "seqrf")) {
    usethis::ui_info("you specified {usethis::ui_code('seqdata')} which are not stored as sequence object\n     and a valid {usethis::ui_code('seqrfobject')}; the {usethis::ui_code('seqdata')};\n     as well as the potentially specified parameters\n     {usethis::ui_code(c('k', 'sortv', 'weighted', 'grp.meth', 'squared', 'pow'))}\n     will be ignored; the plot will be rendered for the {usethis::ui_field('seqrfobject')}")
  }
  if (inherits(seqdata, "stslist") & !is.null(seqrfobject) & 
      !inherits(seqrfobject, "seqrf")) {
    usethis::ui_info("you specified a {usethis::ui_code('seqrfobject')} & {usethis::ui_code('seqdata')};\n    the {usethis::ui_code('seqrfobject')} is not of class {usethis::ui_code('seqrf')} and will be ignored;\n    the plot will be rendered for the {usethis::ui_field('seqdata')} if the other parameters are specified correctly")
  }
  if (!is.null(seqdata) & !inherits(seqdata, "stslist") & !inherits(seqdata, 
                                                                    "seqrf") & !inherits(seqrfobject, "seqrf")) {
    stop("you specified seqdata which are not stored as sequence object and no valid seqrfobject;\nuse 'TraMineR::seqdef' to create a sequence object of class 'stslist' or specify a valid seqrfobject)")
  }
  if (!is.null(seqdata) & !inherits(seqdata, "stslist") & inherits(seqdata, 
                                                                   "seqrf")) {
    stop("you specified seqdata which are of class 'seqrf';\nprobably you forgot to type 'seqrfobject = '")
  }
  if (is.null(seqrfobject) & (is.null(seqdata) | !inherits(seqdata, 
                                                           "stslist"))) {
    stop("no seqrfobject specified & seqdata are either not specified or not\nstored as sequence object; use 'TraMineR::seqdef' to create one")
  }
  if (!inherits(seqrfobject, "seqrf") & is.null(diss)) {
    stop("no seqrfobject specified & diss = NULL; provide a dissimilarity matrix\nprovide a dissimilarity matrix ('diss')")
  }
  if (is.null(border)) 
    border <- FALSE
  if (!is.logical(yaxis) | !is.logical(quality)) {
    stop("the arguments `yaxis`, and `quality`  have to be objects of type logical")
  }
  if (which.plot %in% c("both", "medoids", "diss.to.med") == 
      FALSE) {
    stop("`which.plot` must take one of the following values: \"both\", \"medoids\", \"diss.to.med\"")
  }
  if (!inherits(seqrfobject, "seqrf")) {
    seqrfobject <- TraMineR::seqrf(seqdata = seqdata, diss = diss, 
                                   k = k, sortv = sortv, weights = NULL, weighted = weighted, 
                                   grp.meth = grp.meth, squared = squared, pow = pow)
  }
  seqdata <- seqrfobject$seqtoplot
  k <- nrow(seqdata)
  if (is.null(ylab)) 
    ylab <- "Frequency group"
  ylabels <- pretty(1:k)
  ylabels[1] <- 1
  ylabels[length(ylabels)] <- k
  if (ylabels[length(ylabels)] == ylabels[length(ylabels) - 
                                          1] + 1) {
    ylabels <- ylabels[ylabels != ylabels[length(ylabels) - 
                                            1]]
  }
  if (ylabels[1] == ylabels[2] - 1) {
    ylabels <- ylabels[ylabels != ylabels[2]]
  }
  ybrks <- ylabels
  if (is.null(box.color)) 
    box.color <- "black"
  if (is.null(box.fill)) 
    box.fill <- "white"
  if (is.null(box.alpha)) 
    box.alpha <- 1
  if (is.null(outlier.color)) 
    outlier.color <- "black"
  if (is.null(outlier.fill)) 
    outlier.fill <- "transparent"
  if (is.null(outlier.alpha)) 
    outlier.alpha <- 1
  if (outlier.jitter.height > 0.375) 
    outlier.jitter.height <- 0.375
  suppressMessages(p1 <- ggseqiplot(seqdata, border = border, 
                                    weighted = FALSE) + labs(title = "Group medoids", y = ylab) + 
                     theme(plot.title = element_text(hjust = 0.5)))
  ylabsbrks <- ggplot_build(p1)$layout$panel_params[[1]]$y$get_labels()
  if (which.plot == "medoids") {
    suppressMessages(p1 <- ggseqiplot(seqdata, border = border) + 
                       labs(y = ylab))
  }
  dlist <- seqrfobject$rf$dist.list
  dweights <- seqrfobject$rf$weights.list
  g <- length(dlist)
  wtd.fivenum.tmr <- utils::getFromNamespace("wtd.fivenum.tmr", 
                                             "TraMineR")
  if (inherits(seqrfobject, "seqrfprop")) {
    dist.stat <- matrix(rep(NA, 5 * g), nrow = 5)
    for (i in 1:g) {
      dist.stat[, i] <- wtd.fivenum.tmr(dlist[[i]], dweights[[i]])
    }
  }
  else {
    dist.stat <- sapply(dlist, stats::fivenum)
  }
  rownames(dist.stat) <- c("minimum", "lower-hinge", "median", 
                           "upper-hinge", "maximum")
  boxdata <- dplyr::select(dplyr::mutate(dplyr::mutate(dplyr::as_tibble(t(dist.stat)), 
                                                       k = dplyr::row_number(), .before = 1), aux = 1.5 * (.data$`upper-hinge` - 
                                                                                                             .data$`lower-hinge`), minimum = .data$`lower-hinge` - 
                                           .data$aux, minimum = ifelse(.data$minimum < 0, 0, .data$minimum), 
                                         aux2 = .data$maximum, maximum = .data$`upper-hinge` + 
                                           .data$aux), 1:6)
  colnames(boxdata) <- c("k", "ymin", "lower", "middle", "upper", 
                         "ymax")
  dotdata <- dplyr::filter(dplyr::left_join(dplyr::bind_rows(purrr::imap(seqrfobject$rf$dist.list, 
                                                                         ~dplyr::tibble(k = .y, values = .x))), dplyr::select(boxdata, 
                                                                                                                              .data$k, .data$ymin, .data$ymax), by = "k"), .data$values < 
                             .data$ymin | .data$values > .data$ymax)
  boxdata <- dplyr::relocate(dplyr::left_join(dplyr::summarise(dplyr::group_by(dplyr::mutate(dplyr::left_join(dplyr::bind_rows(purrr::imap(seqrfobject$rf$dist.list, 
                                                                                                                                           ~dplyr::tibble(k = .y, values = .x))), boxdata, by = "k"), 
                                                                                             aux_max = .data$values > .data$upper & .data$values <= 
                                                                                               .data$ymax, aux_min = .data$values < .data$lower & 
                                                                                               .data$values >= .data$ymin, aux_max = ifelse(.data$aux_max == 
                                                                                                                                              TRUE, .data$values, .data$upper), aux_min = ifelse(.data$aux_min == 
                                                                                                                                                                                                   TRUE, .data$values, .data$lower)), .data$k), ymin = min(.data$aux_min), 
                                                               ymax = max(.data$aux_max)), dplyr::select(boxdata, -c(.data$ymin, 
                                                                                                                     .data$ymax)), by = "k"), .data$ymax, .after = .data$upper)
  p2 <- ggplot() + geom_jitter(data = dotdata, aes(y = .data$k, 
                                                   x = .data$values), height = outlier.jitter.height, width = 0, 
                               color = outlier.color, fill = outlier.fill, shape = outlier.shape, 
                               size = outlier.size, stroke = outlier.stroke, alpha = outlier.alpha) + 
    geom_boxplot(data = boxdata, aes(group = .data$k, y = .data$k, 
                                     xmin = .data$ymin, xlower = .data$lower, xmiddle = .data$middle, 
                                     xupper = .data$upper, xmax = .data$ymax), stat = "identity", 
                 orientation = "y", width = 0.75, color = box.color, 
                 fill = box.fill, alpha = box.alpha) + scale_y_continuous(limits = range(0, 
                                                                                         k) + 0.5, breaks = ylabsbrks, expand = expansion(add = c(0, 
                                                                                                                                                  0))) + labs(title = "Dissimilarities to medoid", x = "", 
                                                                                                                                                              y = "") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5), 
                                                                                                                                                                                                axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                                                                                                                                                                                                axis.title.y = element_blank(), axis.line.x = element_line(linewidth = 0.3), 
                                                                                                                                                                                                axis.ticks.x = element_line(linewidth = 0.3), panel.grid.major = element_blank(), 
                                                                                                                                                                                                panel.grid.minor = element_blank())
  if (which.plot == "diss.to.med") {
    p2 <- p2 + scale_y_continuous(limits = range(0, k) + 
                                    0.5, breaks = ylabsbrks, expand = expansion(add = 0)) + 
      labs(y = ylab, x = "", title = NULL) + theme_minimal() + 
      theme(axis.ticks.y = element_blank(), axis.line.x = element_line(linewidth = 0.3), 
            axis.ticks.x = element_line(linewidth = 0.3), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  ggrfplot <- p1 + p2 + plot_layout(widths = c((1 + sqrt(5))/2, 
                                               1), guides = "collect") & theme(legend.position = "bottom")
  if (yaxis == FALSE) {
    ggrfplot[[1]] <- ggrfplot[[1]] + theme(axis.text.y = element_blank(), 
                                           axis.ticks.y = element_blank(), axis.title.y = element_blank())
  }
  if (which.plot == "medoids") 
    ggrfplot <- p1
  if (which.plot == "diss.to.med") 
    ggrfplot <- p2
  if (quality == TRUE) {
    ggrfplot <- ggrfplot + patchwork::plot_annotation(caption = paste("Representation quality: R2 =", 
                                                                      format(round(seqrfobject$rf$R2, 2), nsmall = 2), 
                                                                      "and F =", format(round(seqrfobject$rf$Fstat, 2), 
                                                                                        nsmall = 2)))
  }
  ggrfplot <- ggrfplot + theme(plot.margin = margin(15, 15, 
                                                    10, 15))
  ggrfplot$rfsummary <- summary(seqrfobject)[c(1, 3, 4, 5)]
  return(ggrfplot)
}
