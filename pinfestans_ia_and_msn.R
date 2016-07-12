## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      fig.width = 10,
                      fig.height = 10)
options(width = 100)

## ----packages, message = FALSE-------------------------------------------
library('poppr')
library('ggplot2')
library('RColorBrewer')

## ----data_setup----------------------------------------------------------
pinf <- read.genalex("pinfestans_mx_sa.csv", ploidy = 3)
splitStrata(pinf) <- ~Continent/Country
setPop(pinf) <- ~Continent
pinfreps <- c(Pi02 = 2, D13 = 2, Pi33 = 6, Pi04 = 2, Pi4B = 2, Pi16 = 2,
              G11 = 2, Pi56 = 2, Pi63 = 3, Pi70 = 3, Pi89 = 2)
pinfreps <- fix_replen(pinf, pinfreps)

## ----ia_analysis, fig.show='hide', results = 'asis'----------------------
set.seed(20160730)
res <- poppr(pinf, sample = 999, quiet = TRUE,
             H = FALSE, G = FALSE, lambda = FALSE, E5 = FALSE,
             total = FALSE)
set.seed(20160730)
res_cc <- poppr(pinf, sample = 999, clonecorrect = TRUE, quiet = TRUE,
             H = FALSE, G = FALSE, lambda = FALSE, E5 = FALSE,
             strata = ~Continent/Country, keep = 1, total = FALSE)
knitr::kable(res, digits = 2, caption = "Total Data Set")
knitr::kable(res_cc, digits = 2, caption = "Clone Corrected Data")
p <- last_plot()
op <- p # save the original plot

## ----ia_plot-------------------------------------------------------------
p <- p + facet_wrap(~population, ncol = 1, scales = "free")
p <- p + theme_bw() # initial theme
p <- p + theme(text = element_text(size = rel(5))) # make everything bigger
p <- p + theme(plot.title = element_text(size = rel(4))) # make the title bigger
p <- p + theme(axis.ticks.y = element_blank()) # remove y axis
p <- p + theme(axis.text.y = element_blank())  # remove y axis
p <- p + theme(strip.background = element_blank()) # remove facet label backgrounds
p <- p + theme(strip.text = element_text(face = "bold")) # make facet labels bold
p <- p + theme(panel.grid = element_blank()) # remove gridlines
p <- p + labs(title = "Index of Assocation for clone-corrected data")
p <- p + labs(y = NULL)
p

## ----msn_plot------------------------------------------------------------
setPop(pinf) <- ~Country
min_span_net <- bruvo.msn(pinf, replen = pinfreps, add = TRUE, loss = FALSE,
                          showplot = FALSE, include.ties = TRUE)
set.seed(70)
plot_poppr_msn(pinf,
               min_span_net,
               inds = "none",
               mlg = FALSE,
               gadj = 6,
               nodebase = 1.15,
               palette = RColorBrewer::brewer.pal(4, "Set1"),
               cutoff = NULL,
               quantiles = FALSE,
               beforecut = TRUE)

## ----session_information-------------------------------------------------
devtools::session_info()

