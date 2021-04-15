#!/usr/bin/env Rscript
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly=TRUE)
# check for -h or --help
if ((length(args) == 0) || (any(grep("^(--help|-h)$", args))))
{
    cat("usage: ./plot.R in.csv", sep="\n")
    quit("no", 1)
}

res <- read.csv(args[[1]], header=TRUE, stringsAsFactors=FALSE)
basename <- tools::file_path_sans_ext(basename(args[[1]]))

# calculate the number of states
states <- unique(res$State)
# remove SYM if states are 2 (it's identical to ER)
if (length(states) == 2) {
    res <- res[res$Model != 'SYM', ]
}



pua <- ggplot(res[res$Node == 'PUA', ], aes(x=Probability, group=State, fill=State)) +
    geom_histogram() +
    scale_fill_brewer(palette = 'Set1') +
    facet_grid(State~Model) +
    guides(fill = "none") +
    theme_bw() + ggtitle("PUA")


pnua <- ggplot(res[res$Node == 'PNUA', ], aes(x=Probability, group=State, fill=State)) +
    geom_histogram() +
    scale_fill_brewer(palette = 'Set1') +
    facet_grid(State~Model) +
    guides(fill = "none") +
    theme_bw() + ggtitle("PNUA")


psua <- ggplot(res[res$Node == 'PSUA', ], aes(x=Probability, group=State, fill=State)) +
    geom_histogram() +
    scale_fill_brewer(palette = 'Set1') +
    facet_grid(State~Model) +
    guides(fill = "none") +
    theme_bw() + ggtitle("PSUA")


ggsave(pua, filename=paste0(basename, '-pua.pdf'))
ggsave(psua, filename=paste0(basename, '-psua.pdf'))
ggsave(pnua, filename=paste0(basename, '-pnua.pdf'))


# AICc
p <- ggplot(res, aes(x=AICc, fill=Model)) +
    geom_histogram() +
    scale_fill_brewer(palette = 'Set1') +
    facet_grid(Model~.) +
    guides(fill = "none") +
    theme_bw()

ggsave(p, filename=paste0(basename, '-AICc.pdf'))


# LRT
r.ard <- res[res$Model == 'ARD', c("Tree", "LogLikelihood")]
r.sym <- res[res$Model == 'SYM', c("Tree", "LogLikelihood")]
r.er <- res[res$Model == 'ER', c("Tree", "LogLikelihood")]

r.er <- r.er[!duplicated(r.er), ]
r.sym <- r.sym[!duplicated(r.sym), ]
r.ard <- r.ard[!duplicated(r.ard), ]

r.er <- r.er[ order(r.er$Tree), ]
r.sym <- r.sym[ order(r.sym$Tree), ]
r.ard <- r.ard[ order(r.ard$Tree), ]

cat("\n\nMODEL SELECTION\n")

# https://www.r-phylo.org/wiki/HowTo/Ancestral_State_Reconstruction
# binary
if (length(states) == 2) {
    dof <- 1
    lhr <- 1-pchisq(2*abs(r.er$LogLikelihood - r.ard$LogLikelihood), dof)
    sig <- length(lhr[lhr < 0.05])
    cat(sprintf(
        "\nBinary Data ARD > ER (degrees of freedom = %d): %d/%d\n",
        dof, sig, length(lhr)
    ))

} else if (length(states) == 3) {
    # ternary

    # ARD::
    #     C   I   N
    # C   -   1   2
    # I   3   -   4
    # N   5   6   -
    #
    #
    # SYM:
    #     C   I   N
    # C   -   1   2
    # I   1   -   3
    # N   2   3   -
    #
    # sym vs er
    dof <- 2  # 3 - 1
    lhr <- 1-pchisq(2*abs(r.er$LogLikelihood - r.ard$LogLikelihood), dof)
    sig <- length(lhr[lhr < 0.05])
    cat(sprintf(
        "\nBinary Data SYM > ER (degrees of freedom = %d): %d/%d\n",
        dof, sig, length(lhr)
    ))

    # ard vs er
    dof <- 5  # 6 - 1
    lhr <- 1-pchisq(2*abs(r.er$LogLikelihood - r.ard$LogLikelihood), dof)
    sig <- length(lhr[lhr < 0.05])
    cat(sprintf(
        "\nBinary Data ARD > ER (degrees of freedom = %d): %d/%d\n",
        dof, sig, length(lhr)
    ))

    # sym vs ard
    dof <- 3  # 6 - 3
    lhr <- 1-pchisq(2*abs(r.er$LogLikelihood - r.ard$LogLikelihood), dof)
    sig <- length(lhr[lhr < 0.05])
    cat(sprintf(
        "\nBinary Data ARD > SYM (degrees of freedom = %d): %d/%d\n",
        dof, sig, length(lhr)
    ))

} else {
    quit("Wrong state count!")
}

cat("\n\nAIC:\n")
aic <- c()
for (t in unique(res$Tree)) {
    r <- res[res$Tree == t, c("Model", "AICc")]
    r <- r[!duplicated(r), ]
    aic <- c(aic, r[which(r$AICc == min(r$AICc)), 'Model'])
}

for (m in unique(res$Model)) {
    cat(sprintf("  %s\t%d/%d\n", m, length(aic[aic == m]), length(aic)))
}


