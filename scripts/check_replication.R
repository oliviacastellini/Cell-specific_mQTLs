library(dplyr)
library(data.table)
library(tidyr)

standardise <- function(d, ea="ea", oa="oa", beta="beta", chr="chr", pos="pos") {
    toflip <- d[[ea]] > d[[oa]]
    d[[beta]][toflip] <- d[[beta]][toflip] * -1
    temp <- d[[oa]][toflip]
    d[[oa]][toflip] <- d[[ea]][toflip]
    d[[ea]][toflip] <- temp
    d[["snpid"]] <- paste0(d[[chr]], ":", d[[pos]], "_", toupper(d[[ea]]), "_", toupper(d[[oa]]))
    d
}

repl_obs_exp <- function(b_disc, b_rep, se_disc, se_rep, alpha) {
  p_sign <- pnorm(-abs(b_disc) / se_disc) * pnorm(-abs(b_disc) / se_rep) + ((1 - pnorm(-abs(b_disc) / se_disc)) * (1 - pnorm(-abs(b_disc) / se_rep)))
  p_sig <- pnorm(-abs(b_disc) / se_rep + qnorm(alpha / 2)) + (1 - pnorm(-abs(b_disc) / se_rep - qnorm(alpha / 2)))
  p_rep <- pnorm(abs(b_rep)/se_rep, lower.tail=FALSE)
  res <- tibble::tibble(
    nsnp=length(b_disc),
    metric=c("Sign", "Sign", "P-value", "P-value"),
    datum=c("Expected", "Observed", "Expected", "Observed"),
    value=c(sum(p_sign, na.rm=TRUE), sum(sign(b_disc) == sign(b_rep)), sum(p_sig, na.rm=TRUE), sum(p_rep < alpha, na.rm=TRUE))
  )
  return(list(res=res, variants=dplyr::tibble(sig=p_sig, sign=p_sign)))
}

# Read in clumped
load("16_clumped.rdata")

# Organise data
clumped <- ungroup(clumped)
clumped$snpchr <- gsub("chr", "", clumped$snpchr)
clumped2 <- standardise(clumped, "Allele2", "Allele1", "Effect", "snpchr", "snppos")


##Loop to read and process all 100 chunks 
# Initialize an empty data frame to store the results
merged_data <- data.frame()

# Loop through 100 chunks of neu data
for (i in 1:100) {
  # Read in the current chunk
  chunk_file <- paste0("/home/uy23281/Neu_mQTL/sorted/EPIC_mQTL/results/05/SortedNeu.res.cis.", i, ".RData")
  neu <- fread(chunk_file)

  # Subset the current chunk based on the clumped$cpg values
  neu <- subset(neu, gene %in% clumped$cpg)

  # Separate SNP column into chr and pos
  neu <- neu %>% tidyr::separate(SNP, sep = ":", into = c("chr", "pos"))

  # Separate pos column into pos, a1, and a2
  neu <- neu %>% tidyr::separate(pos, sep = "_", into = c("pos", "a1", "a2"))

  # Standardize neu data
  neu2 <- standardise(neu, "a1", "a2", "beta", "chr", "pos")

  # Merge with clumped2 data
  m <- inner_join(clumped2, neu2, by = c("snpid", "cpg" = "gene"))

  # Calculate se_neu
  m$se_neu <- m$beta / m$`t-stat`

  # Append the current chunk's data to the merged_data
  if (i == 1) {
    merged_data <- m
  } else {
    merged_data <- rbind(merged_data, m)
  }
}
# Look at relationship between effects
summary(lm(Effect ~ beta, merged_data))

# Look at overall replication rate
table(merged_data$`p-value` < 0.05) %>% prop.table

# Look at replication rate (observed vs expected)
repl_obs_exp(merged_data$Effect, merged_data$beta, merged_data$StdErr, merged_data$se_neu, 5e-6)

# plot
pdf("/home/uy23281/Neu_mQTL/sorted/EPIC_mQTL/results/05/neu.pdf")
plot(Effect ~ beta, merged_data)
dev.off()

