#Load Packages
library(phyloseq)
library("ggplot2")
library("RColorBrewer")
library("gplots")
library("Heatplus")
library("vegan")
library("ade4")
library("picante")
theme_set(theme_bw())


**********************************************************************************************
///////////////////////////INITIAL SETUP OF DATA FILES/////////////////////////////////////////
**********************************************************************************************

##Load the files needed
file = "otu_table_NTM.Merged.a.biom"
map = "New.Merged.Map.NTM.m.txt" #leave A1 empty in the mapping file

# Load the abundace table and mapping table 
abundance.table = import_biom(file, taxaPrefix=F)
mapping.table=sample_data(read.table(map, header=T, sep="\t", row.names=1))


lung.physeq=phyloseq(otu_table(abundance.table),tax_table(abundance.table), mapping.table)

#Give a colnames to separate different taxonomic levels
colnames(tax_table(lung.physeq))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU")

# Load the tree file (use the unannotated.tree)
treefile = "97_otus_unannotated.tree"
tree.obj = import_qiime(treefilename = treefile) 


# Now merge the three separate phyloseq objects into a single object
NTM.OTU.Table = merge_phyloseq(lung.physeq, mapping.table, tree.obj)


# Remove taxa with 0 abundance
NTM.OTU.Table = subset_taxa(NTM.OTU.Table, rowSums(otu_table(NTM.OTU.Table)) != 0)


##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
    x/sum(x)
}
NTM.OTU.Rel.Table = transformSampleCounts(NTM.OTU.Table, normalizeSample)

#Remove Tech replicates
NTM.OTU.Rel.Table = subset_samples(NTM.OTU.Rel.Table, Technical_Repeat==0)

#Create phyllum and order tables (do it after normalization and out of the relative table)
Phylum.rel.table = tax_glom(NTM.OTU.Rel.Table, taxrank = "Phylum")
Class.rel.table = tax_glom(NTM.OTU.Rel.Table, taxrank = "Class")
Order.rel.table = tax_glom(NTM.OTU.Rel.Table, taxrank = "Order")
Family.rel.table = tax_glom(NTM.OTU.Rel.Table, taxrank = "Family")
Genus.rel.table = tax_glom(NTM.OTU.Rel.Table, taxrank = "Genus")
OTU.rel.table = tax_glom(NTM.OTU.Rel.Table, taxrank = "OTU")

***********************************************************************************
***********************************************************************************
/////////////////////ANALAYSIS OF ORAL WASH & SPUTUM SAMPLES///////////////////////
***********************************************************************************
***********************************************************************************

#We want to describe the oral and sputum samples from this cohort
#Subset samples of interest: All Sputum and oral samples
Sputum.Oral.Samples.OTU.Rel.Table = subset_samples(NTM.OTU.Rel.Table, Sample_Type_Code %in% c(19,21))
Sputum.Oral.Samples.Genus.rel.table = subset_samples(Genus.rel.table, Baseline_1st==1)
rownames(sample_data(Sputum.Oral.Samples.Genus.rel.table))

Sputum.Samples.Genus.rel.table = subset_samples(Sputum.Oral.Samples.Genus.rel.table, Sample_Type_Code ==21)


/////////////////////////////////////////////////////////////////////////////////////////////////////////
*************************************************HEAT MAP************************************************
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#Heat map with clustering for all Sputum and oral samples
##Here we can select for genera present in >3% relative abundance in 0.5% of the samples (this approach brings in > 70% of the data in almost all samples) WE LIKE THIS APPROACH
Sputum.Oral.Samples.Genus.Rel.wh1 = genefilter_sample(Sputum.Oral.Samples.Genus.rel.table, filterfun_sample(function(x) x > 0.03), A = 0.001 * nsamples(Sputum.Oral.Samples.Genus.rel.table))
Sputum.Oral.Samples.Genus.Rel.table1B = prune_taxa(Sputum.Oral.Samples.Genus.Rel.wh1, Sputum.Oral.Samples.Genus.rel.table)
colnames(sample_data(Sputum.Oral.Samples.Genus.Rel.table1B))
plot_bar(Sputum.Oral.Samples.Genus.Rel.table1B, fill="Genus")

#set data tables  
Sputum.Oral.Samples.GenusData <-otu_table(Sputum.Oral.Samples.Genus.Rel.table1B) #pruned to selected Genuses based on abundance


#create vector to label by Type_Sample
colnames(sample_data(Sputum.Oral.Samples.Genus.Rel.table1B))
sample_data(Sputum.Oral.Samples.Genus.Rel.table1B)$Sample_Type_Code
SampleVector = sample_data(Sputum.Oral.Samples.Genus.Rel.table1B)$Sample_Type_Code #There are 2 different ones

#duplicate to create a color vector and replace value w/ color 
#Colorvector can only replace numbers! 
Colorvector <-SampleVector
Colorvector <- replace(Colorvector, which (Colorvector == "19"), "red")
Colorvector <- replace(Colorvector, which (Colorvector == "21"), "blue")


##Cluster Bray Heatmap
#cluster Genuses(row)
Sputum.Oral.Samples.GenusData.GenusData.Bray.dist <-vegdist(Sputum.Oral.Samples.GenusData, method = "bray")
Sputum.Oral.Samples.GenusData.Genus.Bray.clus <-hclust(Sputum.Oral.Samples.GenusData.GenusData.Bray.dist, "aver")

#cluster samples(Col)
Sputum.Oral.Samples.GenusData.Samples.Bray.dist = distance(Sputum.Oral.Samples.GenusData, method="bray")
Sputum.Oral.Samples.GenusData.Samples.cluster.Bray = hclust(Sputum.Oral.Samples.GenusData.Samples.Bray.dist, "aver")

#Here we are able to change the names for genuses that are labelled as "g__" --> Come back to this
tax_table(Sputum.Oral.Samples.Genus.Rel.table1B)
Sputum.Oral.Samples.Genus.Rel.table1B.New.Names = prune_taxa(tail(names(sort(taxa_sums(Sputum.Oral.Samples.Genus.Rel.table1B))), ntaxa(Sputum.Oral.Samples.Genus.Rel.table1B)), Sputum.Oral.Samples.Genus.Rel.table1B)

# Add a new rank, Strain, with the Genus ids
tax_table(Sputum.Oral.Samples.Genus.Rel.table1B.New.Names) <- cbind(tax_table(Sputum.Oral.Samples.Genus.Rel.table1B.New.Names), Strain=taxa_names(Sputum.Oral.Samples.Genus.Rel.table1B.New.Names))
# Define the ranks you want to include
myranks = c("Class", "Order", "Family", "Genus")
mylabels = apply(tax_table(Sputum.Oral.Samples.Genus.Rel.table1B.New.Names)[, myranks], 1, paste, sep="", collapse="_")

# Add concatenated labels as a new rank after strain
tax_table(Sputum.Oral.Samples.Genus.Rel.table1B.New.Names) <- cbind(tax_table(Sputum.Oral.Samples.Genus.Rel.table1B.New.Names), catglab=mylabels)
tax_table(Sputum.Oral.Samples.Genus.Rel.table1B.New.Names)


#Now Plot Heat map with dendograms
mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))

pdf("Bray_Curtis_Heatmap_All.Sputum.Oral.Samples.pdf", height = 10, width = 20)
heatmap.2(Sputum.Oral.Samples.GenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Sputum.Oral.Samples.GenusData.Genus.Bray.clus),
	Colv = as.dendrogram(Sputum.Oral.Samples.GenusData.Samples.cluster.Bray),
	labRow=tax_table(Sputum.Oral.Samples.Genus.Rel.table1B.New.Names)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(Sputum.Oral.Samples.Genus.Rel.table1B)$Paper_ID,
	cexCol = .8,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	ColSideColors=Colorvector,
	main = "Heatmap of Bray-Curtis Distance",
)
dev.off()


/////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************BETA DIVERSITY*************************************************
/////////////////////////////////////////////////////////////////////////////////////////////////////////
##Calculate distance with Weighted Unifrac 
Sputum.Oral.Samples.wUniF.dist = distance(Sputum.Oral.Samples.OTU.Rel.Table, "wunifrac")

#Estimate the number of axes
Sputum.Oral.Samples.wUniF.pco = dudi.pco(cailliez(Sputum.Oral.Samples.wUniF.dist))

#Plot the PCA
pdf(file="Oral.vs.Sputum.Weighted.UniFrac.pdf", width=6, height=4)
s.class(Sputum.Oral.Samples.wUniF.pco$li, interaction(sample_data(Sputum.Oral.Samples.OTU.Rel.Table)$SampleType), col=c("red", "blue"))
dev.off()

#Test The Statistics
adonis(Sputum.Oral.Samples.wUniF.dist ~ SampleType, data=data.frame(sample_data(Sputum.Oral.Samples.OTU.Rel.Table)))
            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
SampleType   1    0.6028 0.60283  15.045 0.04179  0.001 ***
Residuals  345   13.8235 0.04007         0.95821           
Total      346   14.4263                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

/////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************ALPHA DIVERSITY***********************************************
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#Just compare oral and sputum samples
Oral.Sputum.Shannon_diversity = diversity(otu_table(Sputum.Oral.Samples.OTU.Rel.Table), index = "shannon", MARGIN = 2, base = exp(1))
Oral.Sputum.Simpson_diversity = diversity(otu_table(Sputum.Oral.Samples.OTU.Rel.Table), index = "simpson", MARGIN = 2, base = exp(1))

## TO PLOT IT
boxplot(Oral.Sputum.Shannon_diversity ~ sample_data(Sputum.Oral.Samples.OTU.Rel.Table)$SampleType)

#Test Statistics
wilcox.test(Oral.Sputum.Shannon_diversity ~ sample_data(Sputum.Oral.Samples.OTU.Rel.Table)$SampleType) 
W = 17531, p-value = 0.006062
alternative hypothesis: true location shift is not equal to 0

# To write it
write.table(Oral.Sputum.Shannon_diversity, file="Oral.Sputum.Shannon_diversity.txt", sep="\t") 

## Disease ##
# After modifying table to prep data for graphing, import it back to R
# 1st column: titled "ind"
# 2nd column: titled "values"
shannon.div_oral.sputum <-read.table("Oral.Sputum.Shannon_diversity_for.graph.txt", header = TRUE)

# Prep for scatterplot
xcoord <- rep(0, length(shannon.div_oral.sputum$ind))
xcoord[shannon.div_oral.sputum$ind=="Oral.Rinse"]<- 1
xcoord[shannon.div_oral.sputum$ind=="Sputum"]<- 2


#Explore y axis limits
plot(shannon.div_oral.sputum$values ~ jitter(xcoord, 0.75), pch = 20, xlim=c(0.5, 2.5))
#Then plot

pdf(file="Shannon.Diversity.oral.sputum 1.pdf", height = 8, width = 6)
boxplot(Oral.Sputum.Shannon_diversity ~ sample_data(Sputum.Oral.Samples.OTU.Rel.Table)$SampleType, col=c("red", "blue"), outline=FALSE, boxwex=0.4, xlim=c(0.5, 2.5), ylim=c(0.5, 4.5), names=c("Oral Wash", "Sputum")) # boxplot
par(new=T) # allows overlaying of 2 graphs
plot(shannon.div_oral.sputum$values ~ jitter(xcoord, 0.75), pch = 20, xlim=c(0.5, 2.5), ylim=c(0.5, 4.5), axes=F, xlab="", ylab="")
dev.off()






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


***********************************************************************************
***********************************************************************************
///////////////COMPARING SAMPLES BASED ON NTM CULTURE DATA/////////////////////////
***********************************************************************************
***********************************************************************************



/////////////////////////////////////////////////////////////////////////////////////////////////////////
******************************************BETA DIVERSITY*************************************************
/////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////
** SPUTUM SAMPLES ONLY**
///////////////////////
#Subset just sputum Samples
Sputum.Samples.OTU.Rel.Table = subset_samples(Sputum.Oral.Samples.OTU.Rel.Table, Sample_Type_Code %in% c(21))

#need to eliminate n.a for Current_Mycobacteria
Current_Mycobacteria.Sputum.Samples.OTU.Rel.Table = subset_samples(Sputum.Samples.OTU.Rel.Table, Current_Mycobacteria_Code %in% c(0,1))

#need to eliminate patients on NTM treatment
ntmabx_at_sample.Sputum.Samples.OTU.Rel.Table = subset_samples(Current_Mycobacteria.Sputum.Samples.OTU.Rel.Table, ntmabx_at_sample %in% c(0))

#Weighted Unifrac
ntmabx_at_sample.Sputum.Samples.wUniF.dist = distance(ntmabx_at_sample.Sputum.Samples.OTU.Rel.Table, "wunifrac")

#Estimate the number of axes
ntmabx_at_sample.Sputum.Samples.wUniF.pco = dudi.pco(cailliez(ntmabx_at_sample.Sputum.Samples.wUniF.dist))

#Plot the PCA
pdf(file="Mycobacteria Status No Treatment Sputum.Weighted.UniFrac.pdf", width=4, height=4)
s.class(ntmabx_at_sample.Sputum.Samples.wUniF.pco $li, interaction(sample_data(ntmabx_at_sample.Sputum.Samples.OTU.Rel.Table)$Current_Mycobacteria), col=c("darkgreen", "red"))
dev.off()

#Test The Statistics
adonis(ntmabx_at_sample.Sputum.Samples.wUniF.dist ~ Current_Mycobacteria, data=data.frame(sample_data(ntmabx_at_sample.Sputum.Samples.OTU.Rel.Table)))
                      Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
Current_Mycobacteria   1    0.0843 0.084336  1.9709 0.01657  0.081 .
Residuals            117    5.0064 0.042789         0.98343         
Total                118    5.0907                  1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




//////////////////
** ORAL WASH ONLY**
//////////////////

#Create Oral Wash Table
Oral.Wash.Samples.OTU.Rel.Table = subset_samples(Sputum.Oral.Samples.OTU.Rel.Table, Sample_Type_Code %in% c(19))

#need to eliminate n.a for Current_Mycobacteria
Current_Mycobacteria.Oral.Wash.Samples.OTU.Rel.Table = subset_samples(Oral.Wash.Samples.OTU.Rel.Table, Current_Mycobacteria_Code %in% c(0,1))

#need to eliminate patients on NTM treatment
ntmabx_at_sample.Oral.Wash.Samples.OTU.Rel.Table = subset_samples(Current_Mycobacteria.Oral.Wash.Samples.OTU.Rel.Table, ntmabx_at_sample %in% c(0))

#Weighted Unifrac
ntmabx_at_sample.Oral.Wash.Samples.wUniF.dist = distance(ntmabx_at_sample.Oral.Wash.Samples.OTU.Rel.Table, "wunifrac")

#Estimate the number of axes
ntmabx_at_sample.Oral.Wash.Samples.wUniF.pco = dudi.pco(cailliez(ntmabx_at_sample.Oral.Wash.Samples.wUniF.dist))

#Plot the PCA
pdf(file="Mycobacteria Treatment Oral.Wash.Weighted.UniFrac.pdf", width=4, height=4)
s.class(ntmabx_at_sample.Oral.Wash.Samples.wUniF.pco $li, interaction(sample_data(ntmabx_at_sample.Oral.Wash.Samples.OTU.Rel.Table)$Current_Mycobacteria), col=c("darkgreen", "red"))
dev.off()

#Test the Statistics
adonis(ntmabx_at_sample.Oral.Wash.Samples.wUniF.dist ~ Current_Mycobacteria, data=data.frame(sample_data(ntmabx_at_sample.Oral.Wash.Samples.OTU.Rel.Table)))
                      Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
Current_Mycobacteria   1    0.0716 0.071578  2.3617 0.01826  0.043 *
Residuals            127    3.8490 0.030307         0.98174         
Total                128    3.9206                  1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




////////////////////////////////////////////////////////////////////////////////////
**********************************ALPHA DIVERSITY**********************************
////////////////////////////////////////////////////////////////////////////////////

///////////////////////
**SPUTUM SAMPLES ONLY**
//////////////////////

#testing alpha based on current mycobacterium culture
ntmabx_at_sample.ALL.Sputum.Shannon_diversity = diversity(otu_table(ntmabx_at_sample.Sputum.Samples.OTU.Rel.Table), index = "shannon", MARGIN = 2, base = exp(1))


# To write it to file
write.table(ntmabx_at_sample.ALL.Sputum.Shannon_diversity, file="ntmabx_at_sample.ALL.NoTx.Sputum.Shannon.diversity.txt", sep="\t")
ntmabx_at_sample.ALL.NoTx.shannon.div_sputum <-read.table("ntmabx_at_sample.ALL.NoTx.Sputum.Shannon.diversity_for_graph.txt", header = TRUE)
xcoord <- rep(0, length(ntmabx_at_sample.ALL.NoTx.shannon.div_sputum $ind))
xcoord[ntmabx_at_sample.ALL.NoTx.shannon.div_sputum $ind=="NTM-ve"]<- 1
xcoord[ntmabx_at_sample.ALL.NoTx.shannon.div_sputum $ind=="NTM+ve"]<- 2

#Plot Shannon Diversity
pdf(file="Shannon.Diversity.sputum.NTM.All.NoTx 1.pdf", height = 8, width = 6)
boxplot(ntmabx_at_sample.ALL.Sputum.Shannon_diversity ~ sample_data(ntmabx_at_sample.Sputum.Samples.OTU.Rel.Table)$Current_Mycobacteria, col=c("chartreuse3", "red"), outline=FALSE, boxwex=0.4, xlim=c(0.5, 2.5), ylim=c(0.5, 4.5)) # boxplot
par(new=T) # allows overlaying of 2 graphs
plot(ntmabx_at_sample.ALL.NoTx.shannon.div_sputum $values ~ jitter(xcoord, 0.75), pch = 20, xlim=c(0.5, 2.5), ylim=c(0.5, 4.5), axes=F, xlab="", ylab="")
dev.off()

#Test The Statistics
wilcox.test(ntmabx_at_sample.Sputum.Shannon_diversity ~ sample_data(ntmabx_at_sample.Sputum.Samples.OTU.Rel.Table)$Current_Mycobacteria) 
W = 1456, p-value = 0.04702
alternative hypothesis: true location shift is not equal to 0

///////////////////////
**ORAL WASH ONLY**
//////////////////////

#testing alpha based on current mycobacterium culture
ntmabx_at_sample.nona.NoTx.Oral.Wash.Shannon_diversity = diversity(otu_table(Current_Mycobacteria.NoTx.Oral.Wash.Samples.Baseline.OTU.Rel.Table), index = "shannon", MARGIN = 2, base = exp(1))

# To write it to file
write.table(ntmabx_at_sample.nona.NoTx.Oral.Wash.Shannon_diversity, file="ntmabx_at_sample.nona.NoTx.Oral.Wash.Shannon.diversity.txt", sep="\t")
ntmabx_at_sample.nona.NoTx.shannon.div_oral_wash <-read.table("ntmabx_at_sample.nona.NoTx.Oral.Wash.Shannon.diversity_for_graph.txt", header = TRUE)
xcoord <- rep(0, length(ntmabx_at_sample.nona.NoTx.shannon.div_oral_wash $ind))
xcoord[ntmabx_at_sample.nona.NoTx.shannon.div_oral_wash $ind=="NTM+ve"]<- 1
xcoord[ntmabx_at_sample.nona.NoTx.shannon.div_oral_wash $ind=="NTM-ve"]<- 2

#Plot Shannon Diversity
pdf(file="Shannon.Diversity.oral.wash.NTM.Baseline.NoTx 1.pdf", height = 8, width = 6)
boxplot(ntmabx_at_sample.nona.NoTx.Oral.Wash.Shannon_diversity ~ sample_data(Current_Mycobacteria.NoTx.Oral.Wash.Samples.Baseline.OTU.Rel.Table)$CurrentMyco_NTMtx, col=c("chartreuse3", "red"), outline=FALSE, boxwex=0.4, xlim=c(0.5, 2.5), ylim=c(1.5, 4.5)) # boxplot
par(new=T) # allows overlaying of 2 graphs

#Test The Statistics
wilcox.test(ntmabx_at_sample.nona.NoTx.Oral.Wash.Shannon_diversity ~ sample_data(Current_Mycobacteria.NoTx.Oral.Wash.Samples.Baseline.OTU.Rel.Table)$CurrentMyco_NTMtx) 
W = 760, p-value = 0.6946
alternative hypothesis: true location shift is not equal to 0




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



***********************************************************************************
***********************************************************************************
/////////////////////ANALAYSIS OF BRONCHOSCOPY SAMPLES/////////////////////////////
***********************************************************************************
***********************************************************************************

#Subset samples of interest: All Bronchoscopic samples
Bronch.Done.Samples.OTU.Rel.Table = subset_samples(NTM.OTU.Rel.Table, Bronch_Done %in% c(1))
Bronch.Done.Samples.Genus.rel.table = subset_samples(Genus.rel.table, Bronch_Done ==1)



///////////////////////////////////////////////////////////////////////////////////
***************************HEAT MAP***********************************************
///////////////////////////////////////////////////////////////////////////////////

**********************************
////ALL BRONCH RELATED SAMPLES////
**********************************
#Remove lung tissue samples
Bronch.Done.No.Lung.Tissue.Samples.OTU.Rel.Table = subset_samples(Bronch.Done.Samples.OTU.Rel.Table, Sample_Type_Code !=7)
Bronch.Done.Samples.No.Lung.Tissue.Genus.rel.table = subset_samples(Bronch.Done.Samples.Genus.rel.table, Sample_Type_Code !=7)

#Remove Patients on Treatment
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.OTU.Rel.Table = subset_samples(Bronch.Done.No.Lung.Tissue.Samples.OTU.Rel.Table, CurrentMyco_NTMtx_code %in% c(1,3))
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.rel.table = subset_samples(Bronch.Done.Samples.No.Lung.Tissue.Genus.rel.table, CurrentMyco_NTMtx_code %in% c(1,3))

#Prune data based on Abundance
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.wh1 = genefilter_sample(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.rel.table, filterfun_sample(function(x) x > 0.01), A = 0.05 * nsamples(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.rel.table))
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B = prune_taxa(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.wh1, Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.rel.table)
colnames(sample_data(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B))
plot_bar(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B, fill="Genus")
tax_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B)

#create genusdata table
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData <-otu_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B) #pruned to selected Genuses based on abundance

##Cluster Bray Heatmap
#cluster Genuses(row)
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData.GenusData.Bray.dist <-vegdist(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData, method = "bray")
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData.Genus.Bray.clus <-hclust(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData.GenusData.Bray.dist, "aver")

#cluster samples(Col)
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData.Samples.Bray.dist = distance(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData, method="bray")
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData.Samples.Bray.Bray = hclust(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData.Samples.Bray.dist, "aver")

#Here we are able to change the names for genuses that are labelled as "g__" --> Come back to this
tax_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B)
Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names = prune_taxa(tail(names(sort(taxa_sums(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B))), ntaxa(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B)), Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B)

# Add a new rank, Strain, with the Genus ids
tax_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names) <- cbind(tax_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names), Strain=taxa_names(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names))
# Define the ranks you want to include
myranks = c("Class", "Order", "Family", "Genus")
mylabels = apply(tax_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names)[, myranks], 1, paste, sep="", collapse="_")

# Add concatenated labels as a new rank after strain
tax_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names) <- cbind(tax_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names), catglab=mylabels)
tax_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names)


#create vector to lable by Type_Sample
colnames(sample_data(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names))
sample_data(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names)$Sample_Type_Code
SampleVector = sample_data(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names)$Sample_Type_Code #There are 7 different ones

#duplicate to create a color vector and replace value w/ color 
#Colorvector can only replace numbers! 
Colorvector <-SampleVector
Colorvector <- replace(Colorvector, which (Colorvector == "1"), "grey")
Colorvector <- replace(Colorvector, which (Colorvector == "4"), "orange")
Colorvector <- replace(Colorvector, which (Colorvector == "5"), "firebrick3")
Colorvector <- replace(Colorvector, which (Colorvector == "6"), "darkturquoise")
Colorvector <- replace(Colorvector, which (Colorvector == "19"), "red")
Colorvector <- replace(Colorvector, which (Colorvector == "21"), "blue")


#Choose Color Palatte
mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))

#Heatmap of samples clustered by sample and taxa, but color coded by DMM Clustering
pdf("Bray_Curtis_Heatmap_No.Lung.Tissue.NoTx.withBronch.pdf", height = 10, width = 20)
	heatmap.2(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData.Genus.Bray.clus),
	Colv = as.dendrogram(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.GenusData.Samples.Bray.Bray),
	labRow=tax_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B.New.Names)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.Genus.Rel.table1B)$Paper_ID,
	cexCol = .8,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	ColSideColors=Colorvector,
	main = "Heatmap of Bray-Curtis Distance",
)
dev.off()



///////////////////////////////////////////////////////////////////////////////////
***************************BETA DIVERSITY******************************************
///////////////////////////////////////////////////////////////////////////////////

#Weighted Unifrac
Bronch.Done.F.wUniF.dist = distance(Bronch.Done.No.Lung.Tissue.Samples.OTU.Rel.Table, "wunifrac")

#Examine Axes
Bronch.Done.F.wUniF.pco = dudi.pco(cailliez(Bronch.Done.F.wUniF.dist))

#Plot PCA
pdf(file="Bronch.Samples.Weighted.UniFrac.pdf", width=4, height=4)
s.class(Bronch.Done.F.wUniF.pco$li, interaction(sample_data(Bronch.Done.No.Lung.Tissue.Samples.OTU.Rel.Table)$SampleType), col=c("darkturquoise",  "grey", "orange",  "red", "blue", "firebrick3"))
dev.off()

#Analyze just BAL samples
BAL.Samples.OTU.Rel.Table = subset_samples(Bronch.Done.No.Lung.Tissue.Samples.OTU.Rel.Table, Sample_Type_Code ==6)

#Remove BAL sample with no culture available
BAL.Samples.NTMCult.Avail.OTU.Rel.Table = subset_samples(BAL.Samples.OTU.Rel.Table, NTM_culture_avail !=0)

#Remove Patients on Treatment
BAL.Samples.NTMCult.NoTx.Avail.OTU.Rel.Table = subset_samples(BAL.Samples.NTMCult.Avail.OTU.Rel.Table, CurrentMyco_NTMtx_code %in% c(1,3))

#Weighted Unifrac
BAL.NTMCult.Avail.NoTx.wUniF.dist = distance(BAL.Samples.NTMCult.NoTx.Avail.OTU.Rel.Table, "wunifrac")

#Examine Axes
BAL.NTMCult.Avail.NoTx.wUniF.pco = dudi.pco(cailliez(BAL.NTMCult.Avail.NoTx.wUniF.dist))

#Plot PCA
pdf(file="BAL.culture.avail.NoTx.Weighted.UniFrac.pdf", width=4, height=4)
s.class(BAL.NTMCult.Avail.NoTx.wUniF.pco$li, interaction(sample_data(BAL.Samples.NTMCult.NoTx.Avail.OTU.Rel.Table)$Current_Mycobacteria), col=c("darkgreen", "red"))
dev.off()

#Test Statistics
adonis(BAL.NTMCult.Avail.NoTx.wUniF.dist ~ Current_Mycobacteria, data=data.frame(sample_data(BAL.Samples.NTMCult.NoTx.Avail.OTU.Rel.Table)))
                     Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
Current_Mycobacteria  1   0.17633 0.176329  2.0465 0.05678  0.051 .
Residuals            34   2.92942 0.086159         0.94322         
Total                35   3.10575                  1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


///////////////////////////////////////////////////////////////////////////////////
*********************************ALPHA DIVERSITY**********************************
///////////////////////////////////////////////////////////////////////////////////

#Create diversity table for Bronch Done Patients
Bronch.Done.No.Lung.Tissue.Shannon_diversity = diversity(otu_table(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.OTU.Rel.Table), index = "shannon", MARGIN = 2, base = exp(1))

# To write it
write.table(Bronch.Done.No.Lung.Tissue.Shannon_diversity, file="Bronch.Done.No.Lung.Tissue.Shannon.diversity.txt", sep="\t")


## Disease ##
# After modifying table to prep data for graphing, import it back to R
# 1st column: titled "ind"
# 2nd column: titled "values"
shannon.div_bronch.done.no.lung.tissue <-read.table("Bronch.Done.No.Lung.Tissue.Shannon_diversity_for.graph.txt", header = TRUE)

# Prep for scatterplot
xcoord <- rep(0, length(shannon.div_bronch.done.no.lung.tissue $ind))
xcoord[shannon.div_bronch.done.no.lung.tissue $ind=="BAL"]<- 6 #darkturquoise
xcoord[shannon.div_bronch.done.no.lung.tissue $ind=="BKG"]<- 1 #grey
xcoord[shannon.div_bronch.done.no.lung.tissue $ind=="Nasal.Swab"]<- 2 #orange
xcoord[shannon.div_bronch.done.no.lung.tissue $ind=="Oral.Rinse"]<- 3 #red
xcoord[shannon.div_bronch.done.no.lung.tissue $ind=="Sputum"]<- 4 #blue
xcoord[shannon.div_bronch.done.no.lung.tissue $ind=="Supraglottic"]<- 5 #firebrick3

#reorder labels
sample_data(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.OTU.Rel.Table)$cut<-factor(sample_data(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.OTU.Rel.Table)$SampleType, levels=c("BKG", "Nasal.Swab", "Oral.Rinse", "Sputum" , "Supraglottic", "BAL"), labels=c("BKG", "Nasal\nSwab", "Oral\nWash", "Sputum" , "Supraglottic", "BAL"))


#Plot Shannon Diversity
pdf(file="Shannon.Diversity.Bronch.Done.No.Lung.Tissue.No.Tx 1.pdf", height = 8, width = 6)
boxplot(Bronch.Done.No.Lung.Tissue.Shannon_diversity ~ sample_data(Bronch.Done.No.Lung.Tissue.Samples.No.Tx.OTU.Rel.Table)$cut, at = c(1,2,3,4,5,6), names =c ("BKG", "Nasal\nSwab", "Oral\nWash", "Sputum" , "SUP", "BAL"), cex.axis=0.7, col=c("grey", "orange",  "red", "blue", "firebrick3", "darkturquoise"), outline=FALSE, xlim=c(0.5, 6.5), ylim=c(0.5, 5.5)) # boxplot
par(new=T) # allows overlaying of 2 graphs
plot(shannon.div_bronch.done.no.lung.tissue $values ~ jitter(xcoord, 0.75), pch = 20, xlim=c(0.5, 6.5), ylim=c(0.5, 5.5), axes=F, xlab="", ylab="")
dev.off()

#Create diversity table for BAL Samples
BAL.Samples.NTMCult.Avail.Shannon_diversity = diversity(otu_table(BAL.Samples.NTMCult.Avail.OTU.Rel.Table), index = "shannon", MARGIN = 2, base = exp(1))

# To write it
write.table(BAL.Samples.NTMCult.Avail.Shannon_diversity, file="BAL.Samples.NTMCult.Avail.Shannon.diversity.txt", sep="\t")

#Test Statistics
wilcox.test(BAL.Samples.NTMCult.Avail.Shannon_diversity ~ sample_data(BAL.Samples.NTMCult.Avail.OTU.Rel.Table)$Current_Mycobacteria) 
W = 134, p-value = 0.2228
alternative hypothesis: true location shift is not equal to 0

## Disease ##
# After modifying table to prep data for graphing, import it back to R
# 1st column: titled "ind"
# 2nd column: titled "values"
shannon.div_BAL.Samples.NTMCult.Avail <-read.table("BAL.Samples.NTMCult.Avail.Shannon_diversity_for.graph.txt", header = TRUE)

# Prep for scatterplot
xcoord <- rep(0, length(shannon.div_BAL.Samples.NTMCult.Avail $ind))
xcoord[shannon.div_BAL.Samples.NTMCult.Avail $ind=="NTM-"]<- 1 #darkturquoise
xcoord[shannon.div_BAL.Samples.NTMCult.Avail $ind=="NTM+"]<- 2 #grey


#Then plot
pdf(file="Shannon.Diversity.BAL.Samples.NTMCult.Avail 1.pdf", height = 8, width = 6)
boxplot(BAL.Samples.NTMCult.Avail.Shannon_diversity ~ sample_data(BAL.Samples.NTMCult.Avail.OTU.Rel.Table)$Current_Mycobacteria, names =c ("NTM-ve", "NTM+ve") , col=c("chartreuse3", "red"), outline=FALSE, boxwex=0.4, xlim=c(0.5, 2.5), ylim=c(0.5, 5.5)) # boxplot
par(new=T) # allows overlaying of 2 graphs
plot(shannon.div_BAL.Samples.NTMCult.Avail$values ~ jitter(xcoord, 0.75), pch = 20, xlim=c(0.5, 2.5), ylim=c(0.5, 5.5), axes=F, xlab="", ylab="")
dev.off()