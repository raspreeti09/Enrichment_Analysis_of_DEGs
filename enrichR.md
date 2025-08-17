# Install and load necessary packages
install.packages("enrichR")  # Install EnrichR if not already installed
install.packages("ggplot2")  # For visualization
library(enrichR)
library(ggplot2)

# List available databases in EnrichR
dbs <- listEnrichrDbs()
print(head(dbs))  # View available databases

# Define databases for enrichment analysis
dbs <- c("GO_Biological_Process_2021", 
         "GO_Cellular_Component_2021", 
         "GO_Molecular_Function_2021", 
         "KEGG_2021_Human")  

# Define the gene list
genes <- c("PLA2G2A", "TLR2", "NLRP3", "PRKCD", "PTGS2", "PDGFRA", "IL6ST", "HPGD", "MIF", "CCR2", "MAPK9", "ALPL", "CTSK")  
#Put your DEGs (genes)
# Perform enrichment analysis
enriched_results <- enrichr(genes, dbs)

# View results
head(enriched_results[["GO_Biological_Process_2021"]])  # BP terms
head(enriched_results[["GO_Cellular_Component_2021"]])  # CC terms
head(enriched_results[["GO_Molecular_Function_2021"]])  # MF terms
head(enriched_results[["KEGG_2021_Human"]])             # KEGG pathways

# Save results to CSV
write.csv(enriched_results[["GO_Biological_Process_2021"]], "GO_BP_EnrichR.csv", row.names=FALSE)
write.csv(enriched_results[["GO_Cellular_Component_2021"]], "GO_CC_EnrichR.csv", row.names=FALSE)
write.csv(enriched_results[["GO_Molecular_Function_2021"]], "GO_MF_EnrichR.csv", row.names=FALSE)
write.csv(enriched_results[["KEGG_2021_Human"]], "KEGG_EnrichR.csv", row.names=FALSE)

# Visualization: Barplot for KEGG enrichment
kegg_df <- enriched_results[["KEGG_2021_Human"]]
ggplot(kegg_df, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "KEGG Pathway Enrichment (EnrichR)", x = "Pathway", y = "-log10(Adj P-value)")

# Visualization: Barplot for GO BP enrichment
bp_df <- enriched_results[["GO_Biological_Process_2021"]]
ggplot(bp_df, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  labs(title = "GO Biological Process Enrichment (EnrichR)", x = "Biological Process", y = "-log10(Adj P-value)")

# Optional: Ridgeplot alternative (not directly available in EnrichR)


##Select top 30 pathways only


library(ggplot2)
library(dplyr)  # Load dplyr for data manipulation

# Select top 30 KEGG pathways based on Adjusted P-value
kegg_df <- enriched_results[["KEGG_2021_Human"]] %>%
  arrange(Adjusted.P.value) %>%  # Sort by p-value (ascending)
  slice_head(n = 20)  # Take top 30

ggplot(kegg_df, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "Top 20 KEGG Pathway Enrichment (EnrichR)", x = "Pathway", y = "-log10(Adj P-value)") +
  theme_minimal()


# Select top 20 GO BP terms
bp_df <- enriched_results[["GO_Biological_Process_2021"]] %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 20)

ggplot(bp_df, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  labs(title = "Top 20 GO Biological Process Enrichment (EnrichR)", x = "Biological Process", y = "-log10(Adj P-value)") +
  theme_minimal()

# Select top 20 GO CC terms
bp_df <- enriched_results[["GO_Cellular_Component_2021"]] %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 20)

ggplot(bp_df, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "orange") +
  coord_flip() +
  labs(title = "Top 20 GO Cellular Component Enrichment (EnrichR)", x = "Cellular Component", y = "-log10(Adj P-value)") +
  theme_minimal()

# Select top 20 GO MF terms
bp_df <- enriched_results[["GO_Molecular_Function_2021"]] %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 20)

ggplot(bp_df, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(title = "Top 20 GO Molecular Function Enrichment (EnrichR)", x = "Molecular Function", y = "-log10(Adj P-value)") +
  theme_minimal()


library(stringr)

# Add KEGG pathway IDs (extract from the "Term" column)
kegg_df <- enriched_results[["KEGG_2021_Human"]]

# Check the structure of the Term column
head(kegg_df$Term)

# Extract hsa ID (if present) from the Term column
kegg_df$Pathway_ID <- str_extract(kegg_df$Term, "hsa\\d{5}")

# Save KEGG results to CSV with hsa IDs
write.csv(kegg_df, "KEGG_EnrichR_with_hsa.csv", row.names=FALSE)

