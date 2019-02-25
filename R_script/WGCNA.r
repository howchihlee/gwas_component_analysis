args = commandArgs(trailingOnly=TRUE)

save_path = args[1]
file_to_read = args[2]

data = read.csv(file_to_read)

library('WGCNA')


scale_free_fit_cut_off = 0.8;
soft_power_search_range = c(2:10) 
minModuleSize = 5;
MEDissThres = 0.2


expr = as(data,"matrix")

sft = pickSoftThreshold(expr, powerVector = soft_power_search_range, verbose = 5, 
                        networkType ="signed hybrid", corFnc= "bicor", corOptions = list(maxPOutliers =0.1))
sft_r2 = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]

softpower = min(which(sft_r2 > scale_free_fit_cut_off))
if (is.infinite(softpower)){
  softpower = soft_power_search_range [which.max(sft_r2)]
} else {
  softpower = soft_power_search_range [softpower]
}
# Plot the results:
sizeGrWindow(9, 5)
pdf(file = paste(save_path, "scale_free_fit.pdf", sep = ''), wi = 9, he = 6)

par(mfrow = c(1,2));

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft_r2,
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], sft_r2,
     labels=soft_power_search_range,
     cex=scale_free_fit_cut_off, col="red");
abline(h=scale_free_fit_cut_off, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels=soft_power_search_range, 
     cex=scale_free_fit_cut_off, col="red")
dev.off()


adjacency = adjacency(expr, power = softpower, corFnc="bicor",
                      type="signed hybrid", corOptions = "use = 'p', maxPOutliers = 0.1")
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)



# Calculate eigengenes
MEList = moduleEigengenes(expr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");

plot_dendrogram <- function(METree, save_path){sizeGrWindow(7, 6)
  pdf(file = paste(save_path, "module_tree_not_merged.pdf", sep = ''), wi = 9, he = 6)
  plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  
  dev.off()}

try(plot_dendrogram(METree, save_path))



merge = mergeCloseModules(expr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
pdf(file = paste(save_path, "gene_tree_dendrogram.pdf", sep = ''), wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#write.table(expr, file = sprintf('%s%s', save_path, 'used_expression.csv'), 
#            quote = F, row.names = rownames(expr), col.names =colnames(expr), sep = ',')
#write.table(TOM, file = sprintf('%s%s', save_path, 'TOM.csv'), quote = F, 
#            row.names = rownames(adjacency), col.names = colnames(adjacency), sep = ',')


write.table(dynamicColors, file = sprintf('%s%s', save_path, 'gene_module.csv'), 
            quote = F, row.names = colnames(expr), col.names =F, sep = ',')

write.table(mergedColors, file = sprintf('%s%s', save_path, 'merged_gene_module.csv'), 
            quote = F, row.names = colnames(expr), col.names =F, sep = ',')

write.table(MEs, file = sprintf('%s%s', save_path, 'eigengenes.csv'), quote = F, 
            row.names = F, col.names = colnames(MEs), sep = ',')

write.table(mergedMEs, file = sprintf('%s%s', save_path, 'merged_eigengenes.csv'), quote = F, 
            row.names = F, col.names = colnames(mergedMEs), sep = ',')
