#!/usr/bin/env Rscript
library(ape)
library(phytools)
library(corHMM)

NONUTOAZ <- c("San_Juan_Pueblo_Tewa-36", "Kiowa-35")

NUA <- c(
    "Hopi-16",
    "Tubatulabel-10",
    "Luiseno-14",
    "Cupeno-13",
    "Cahuilla-12",
    "Gabrielino-15",
    "Serrano-11",
    "Comanche-5",
    "Pannamint-3",
    "Shoshoni_Gosiute_Dialect-4",
    "Chemehuevi-7",
    "Southern_Paiute-8",
    "Southern_Ute-9",
    "Kawaiisu-6",
    "Mono-1"
)


SUA <- c(
    "Yaqui-26",
    "Classical_Aztec-29",
    "Huichol-28",
    "Opata-20",
    "Tarahumara-23",
    "Northern_Tepehuan-19"
)


args <- commandArgs(trailingOnly=TRUE)
# check for -h or --help
if ((length(args) == 0) || (any(grep("^(--help|-h)$", args))))
{
    cat("usage: ./get_rates.R trees data results", sep="\n")
    quit("no", 1)
}


trees <- read.nexus(args[[1]])
data <- read.delim(args[[2]], header = FALSE, stringsAsFactors = FALSE)

results <- data.frame(
    Tree = NULL,
    Model = NULL,
    Node = NULL,
    LogLikelihood = NULL,
    AICc = NULL,
    State = NULL,
    Probability = NULL
)



for (tree_idx in 1:length(trees)) {
    tree <- trees[[tree_idx]]
    clades <- list(
        'PUA' = tree$tip.label[tree$tip.label != NONUTOAZ],
        'PNUA' = tree$tip.label[tree$tip.label %in% NUA],
        'PSUA' = tree$tip.label[tree$tip.label %in% SUA]
    )

    for (model in c("ER", "SYM", "ARD")) {
        cat(paste(tree_idx, model), sep="\n")
        r <- corHMM::rayDISC(
            tree, data, ntraits=1, charnum=1, model=model,
            lewis.asc.bias = TRUE,
            node.states = "marginal",  # joint, marginal, scaled
            # (DEFAULT) assumes equal weighting among all possible root states
            #root.p = NULL
            # follows Maddison et al. (2007) and FitzJohn et al. (2009).
            root.p = "maddfitz"
        )

        for (node in names(clades)) {
            n <- getMRCA(tree, clades[[node]]) - Ntip(tree)
            for (s in colnames(r$states)) {
                results <- rbind(results, data.frame(
                    Tree = tree_idx,
                    Model = model,
                    Node = node,
                    LogLikelihood = r$loglik,
                    AICc = r$AICc,
                    State = s,
                    Probability = r$states[n, s]
                ))
            }
        }

    }
}

write.csv(results, args[[3]], quote = FALSE, row.names = FALSE)

