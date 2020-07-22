#Materials and Methods 
We collected and analyzed the several RNA-seq studies of the normal tissue of cornea (n=19 samples), retina (n=310 samples), retinal pigment epithelium (RPE) (n=207 samples) and lung (n = 546 samples) from the NCBI GEO database (GSE77938 and GSE115828) and GTEx database. Individual data sets underwent stringent quality control and normalization, outliers were defined as samples with standardized sample network connectivity Z scores < -2 and were removed. To specifically characterize the biological pathways involved, we performed k-means analysis, identifying several co-expression clusters. The cluster eigengene was defined as the first principal component summarizing the expression patterns of all genes into a single expression profile within a given cluster. The cluster membership for each virus entry factors and disease susceptibility gene was determined. The gene content of function related module (cluster) was characterized using the Metascape gene enrichment analysis tool. Statistical analysis was done using the R project for statistical computing (http://www.r-project.org). 