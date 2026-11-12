# Builds example/gene_bins.csv: 2 toy genes x 30 epigenomes, three chromatin
# patterns per gene (active, poised/bivalent, repressed) with 10% noisy bins.
set.seed(1)
pat <- list(active = c(rep("1_TssA", 3), rep("5_Tx", 12)),
            bivalent = c(rep("14_TssBiv", 3), rep("17_ReprPCWk", 12)),
            repressed = c(rep("16_ReprPC", 3), rep("18_Quies", 12)))
states <- unique(unlist(pat))
rows <- list()
for (g in c("GENE_A", "GENE_B")) for (s in 1:30) {
  st <- pat[[(s - 1) %/% 10 + 1]]
  noisy <- runif(length(st)) < 0.1; st[noisy] <- sample(states, sum(noisy), TRUE)
  rows[[length(rows) + 1]] <- data.frame(gene_id = g, sample_id = sprintf("sample%02d", s),
                                         bin_id = seq_along(st), state = st, bin_weight = 200)
}
write.csv(do.call(rbind, rows), "example/gene_bins.csv", row.names = FALSE)
