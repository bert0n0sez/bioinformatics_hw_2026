#BiocManager::install("limma")

library(Biobase)
library(limma)

geo_id <- 'GSE63885'
exp = read.csv(file = 'GSE63885/expression_for_limma2.csv', header = TRUE, row.names = 'Gene.Symbol')
ann = read.csv(file = 'GSE63885/annotation_for_limma.csv', header = TRUE, row.names = 'X')


exp_set <- ExpressionSet(assayData = as.matrix(exp), phenoData = AnnotatedDataFrame(ann))

slope <- factor(ann$clinical.status.post.1st.line.chemotherapy..cr...complete.response..pr...partial.response..sd...stable.disease..p...progression..ch1, levels = c("pCR",
                                                                                                                                                                     "pNC"), labels = c(1, 0))

pCR=as.integer(as.vector(slope))

iterse <- rep(1, length(slope))

design <- cbind(npCR=iterse,pCR=pCR)

fit <- lmFit(exp_set, design)
fit <- eBayes(fit)
top <-topTable(fit, coef="pCR", adjust="BH", n = Inf)

pval_threshold <- 0.05

top$significant_logFC1 <- ifelse(
  top$P.Val < pval_threshold & abs(top$logFC) > 1,
  "Significant",
  "Not Significant"
)

top$significant_logFC2 <- ifelse(
  top$P.Val < pval_threshold & abs(top$logFC) > 2,
  "Significant",
  "Not Significant"
)

top$significant_logFC3 <- ifelse(
  top$P.Val < pval_threshold & abs(top$logFC) > 3,
  "Significant",
  "Not Significant"
)

sig1 <- sum(top$significant_logFC1 == "Significant")
sig2 <- sum(top$significant_logFC2 == "Significant")
sig3 <- sum(top$significant_logFC3 == "Significant")

cat("Significant genes, logFC > 1:", sig1, "\n")
cat("Significant genes, logFC > 2:", sig2, "\n")
cat("Significant genes, logFC > 3:", sig3, "\n")

write.csv(top, "LIMMA_results.csv")

library(ggplot2)

ggplot(top, aes(x = logFC, y = -log10(adj.P.Val),
                color = significant_logFC3)) +
  
  geom_point(alpha = 0.7, size = 1.5) +
  
  scale_color_manual(values = c("grey", "red")) +
  
  geom_vline(xintercept = c(-3, 3),
             linetype = "dashed",
             color = "blue") +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "green") +
  
  theme_minimal() +
  
  labs(
    title = "Volcano plot (logFC > 3)",
    x = "logFC",
    y = "-log10(P.Val)"
  )

ggplot(top, aes(x = logFC, y = -log10(adj.P.Val),
                color = significant_logFC2)) +
  
  geom_point(alpha = 0.7, size = 1.5) +
  
  scale_color_manual(values = c("grey", "red")) +
  
  geom_vline(xintercept = c(-2, 2),
             linetype = "dashed",
             color = "blue") +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "green") +
  
  theme_minimal() +
  
  labs(
    title = "Volcano plot (logFC > 2)",
    x = "logFC",
    y = "-log10(P.Val)"
  )

ggplot(top, aes(x = logFC, y = -log10(adj.P.Val),
                color = significant_logFC1)) +
  
  geom_point(alpha = 0.7, size = 1.5) +
  
  scale_color_manual(values = c("grey", "red")) +
  
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "blue") +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "green") +
  
  theme_minimal() +
  
  labs(
    title = "Volcano plot (logFC > 1)",
    x = "logFC",
    y = "-log10(P.Val)"
  )