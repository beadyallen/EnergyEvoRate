
#Load up the required packages
library(phyloseq)
library(ggplot2)
library(vegan)
library(reshape2)
library(stringr)
library(ape)

library(dplyr)
library(data.table)

sessionInfo()

#Assuming the data is conatined in files with names matching *.biom.json 
biomlist <- list.files('Inputs/biom_files/', '*.json', full.names = T)
prj_names <- gsub('\\.biom\\.json','\\1',basename(biomlist))
nmgsdirs <- paste('./NMGS_processing_3/',prj_names, sep="")

#Only required if performing NMGS theta calculations
for (p in nmgsdirs)
#Create the base NMGS directory
    dir.create(p, recursive = T)


classify_data <- read.csv('Inputs/Taxa_Classifications_20200818_only_genus.csv', header = T, comment.char = '#', stringsAsFactors = F )
metabolic_energies <- read.csv('Inputs/Energies_with_RET.csv')

met_paths <- unique(classify_data$Metabolism)
names(met_paths) <- met_paths

#Build up Taxa search expressions for them
Taxa_Search_Expressions <- lapply(met_paths, FUN = function(n){
  X = classify_data[classify_data$Metabolism==n,]
  d <- apply(X, 1, function(X) paste(X[[2]],'=="',X[[1]],'"',sep=''))
  searchstr  = paste(unlist(d),collapse = ' | ')
  parse(text=searchstr)
})


PSobjs_full <- sapply(biomlist, function(x) import_biom(x, parseFunction = parse_taxonomy_greengenes))

PSobjs_full

tt <- lapply(PSobjs_full, function(x) {unique(as.vector(tax_table(x)[,"Genus"]))})

write.table(unique(as.vector(unlist(tt))), "Detected_genera.tsv",sep="\t")

#Named list for all results
all_filtered_list = list()

for (i in seq(length(PSobjs_full))){
    curbiom = PSobjs_full[[i]]
    filtered_list <- lapply(Taxa_Search_Expressions, function(x){
                                return(tryCatch(do.call(subset_taxa, list(curbiom,x)),error=function(e) NULL))
                            })
    filtered_list = filtered_list[!sapply(filtered_list, is.null)]
    all_filtered_list[[prj_names[[i]]]] <- filtered_list
                        
}

allprojects <- list()
genusprojects <- list()
for (i in seq(length(all_filtered_list))){
    raw = PSobjs_full[[i]]
    genus_raw  = tax_table(raw)[!is.na(tax_table(raw)[,"Genus"]),]
    flt = all_filtered_list[[i]]
    allprojects[[names(all_filtered_list)[[i]]]] = sum(sapply(flt, function(x) sum(otu_table(x)))) / sum(otu_table(raw)) * 100
    genusprojects[[names(all_filtered_list)[[i]]]] = sum(sapply(flt, function(x) sum(otu_table(x)))) / sum(otu_table(raw)[rownames(genus_raw),]) * 100    
}
t(rbind(as.data.frame(allprojects), as.data.frame(genusprojects)))

#Generate NMGS input files

for (i in seq(length(all_filtered_list))){
    curbiom <- all_filtered_list[[i]]
    for (x in names(curbiom)){
        outfname = paste(nmgsdirs[[i]],'/',x,'.csv', sep = '')
        write.csv(curbiom[[x]]@otu_table, outfname)
    }
}

#Set up the plot axes - for log scales
xbreaks = c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100),seq(1000,10000,1000))
xlabs <- c("0.001",rep("",9),"0.01",rep("",9),"0.1",rep("",9),"1",rep("",9), "10", rep("",9), "100", rep("",9), "1000", rep("",9))

ybreaks = c(seq(0.1,1,0.1),seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))
ylabs <- c("0.1",rep("",9),"1",rep("",9),"10",rep("",9),"100",rep("",9),"1000",rep("",9))

#Metabolic colours
metab_color_list <-  c("wheat3", "blue","red", "chartreuse3",
                       "darkgoldenrod2","orange","darkred", "purple",
                       "darkgreen", "yellow4", "steelblue3","firebrick2",
                       "indianred", "darkolivegreen4", "honeydew4","cyan", "black")
names(metab_color_list) <- metabolic_energies$Metabolism
metab_label_list = as.character(metabolic_energies$Name)
names(metab_label_list) <- names(metab_color_list)


#Per study labels and colors
study_info_df <- data.frame(
    study=c('AmericanGut_ERP012803', 'Atlantic_ERP012887', 'CanadaWater_ERP012927', 'Donna_ADSludge', 'DuckWater_ERP012631', 'Hydrothermal_ERP011826', 'Lise_Soils', 'MalawianChildren_ERP005437', 'Marine_ERP013553', 'MexicanSoil_SRP037963', 'ParkGrass_SRP044877', 'SevernTrent', 'TARA_ERP001736', 'Volcanic_ERP010094'),
    biome = c("Human","Marine","Freshwater","WWT","Freshwater","Volcanic","Soil", "Human","Marine","Soil","Soil","WWT","Marine", "Volcanic"),
    studylabel = c("American Gut","Atlantic", "Canadian lakes", "Anaerobic sludge", "Farm water", "Hydrothermal vent", "Arctic soil", "Malawian children", "Baltic", "Mexican soil","Canadian parks", "Severn Trent WWT", "TARA Ocean Survey", "Italian Volcanoes")
)

study_name_list = as.character(study_info_df$studylabel)
names(study_name_list) <- as.character(study_info_df$study)
study_biome_list = as.character(study_info_df$biome)
names(study_biome_list) <- as.character(study_info_df$study)
biome_color_list <- c("salmon3","darkseagreen","cornflowerblue","bisque3","darkorange","darkslategray")
names(biome_color_list) <- c("Human", "Marine", "Freshwater", "WWT", "Volcanic", "Soil")


#Carbon source colors + labels
c_src_colours = c("blue","orange")
names(c_src_colours) = c("C","S")

c_src_labels = c("Complex","Simple")
names(c_src_colours) = c("C","S")


#Load the tree
gg_tree = read_tree_greengenes('Inputs/gg_13_5_otus_99_annotated.tree')
knowntips = gg_tree$tip.label

#load the database entries
gg_dat = read.table('Inputs/gg_13_5_taxonomy.txt', sep="\t", header=F, stringsAsFactors = F)
colnames(gg_dat) <- c("taxid","taxstr")
tree_dat = gg_dat[gg_dat$taxid %in% knowntips,]

m = match(gg_tree$tip.label, rownames(tree_dat))
prunetips = gg_tree$tip.label[is.na(head(m))]
length(prunetips)

taxdf = as.data.frame(t(as.data.frame(strsplit(gsub('.__', '', as.character(tree_dat$taxstr)), ';'), stringsAsFactors=FALSE)), stringsAsFactors=FALSE)
colnames(taxdf) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxdf$taxid <- tree_dat$taxid
rownames(taxdf) <- taxdf$taxid

cleantaxfs = sapply(taxdf,function(x) gsub(' ','',x))

cleandf = data.frame(cleantaxfs, stringsAsFactors = F)
rownames(cleandf) = rownames(taxdf)

#Get lists of tip labels we need to search for
tax_groups = lapply(Taxa_Search_Expressions, function(x) rownames(subset(cleandf, eval(x))))

sapply(tax_groups, length)
sum(sapply(tax_groups, length))

ul = unique(unlist(tax_groups))
m = match(gg_tree$tip.label, ul, )
prunelist = gg_tree$tip.label[is.na(m)]

nvec=names(tax_groups)
names(nvec)=nvec
all_prune_taxa = lapply(nvec, function(x) unique(unlist(tax_groups[-which(nvec == x)])))

working_tree = drop.tip(gg_tree, prunelist)

all_trees = lapply(all_prune_taxa, function(x) drop.tip(working_tree,x))

pd_list = lapply(all_trees, function(x) sum(x$edge.length))
pd_df = as.data.frame(t(data.frame(pd_list)))
pd_df[,2] = rownames(pd_df)
colnames(pd_df) <- c("pd","metab")

fulldf = merge(metabolic_energies, pd_df, by.x="Metabolism", by.y="metab")
pd_fulldf <- fulldf[order(fulldf$yield, decreasing = T),]

pd_fulldf 

pd_fits <- pd_fulldf %>% filter(pd > 0.0) %>% do(data.frame(
        yield=summary(lm(log(pd) ~ log(yield), data=.))$coefficients,
        dgcat=summary(lm(log(pd) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

pd_fits$param = c("intercept","slope")

pd_fitsdf<- melt(pd_fits, id.vars = c("param")) %>% dcast(. ~ param +variable)

pd_fulldf %>% filter(pd > 0.0) -> pd2

yield=lm(log(pd) ~ log(yield), data=pd2)
dgcat=lm(log(pd) ~ log(-dgcat), data=pd2)

summary(yield)

summary(dgcat)

pd_plot_yield <- ggplot(pd_fulldf, aes(x=yield,y=pd)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=pd_fitsdf, aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate), 
                  alpha=1.0, color="black",size=2) +
      
      geom_point(size=5, aes(color=Metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      #scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +
      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +
      #scale_color_manual(values = c_src_colours, labels = c_src_labels, name="Anabolic C source" ) +
      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab("Phylogenetic diversity (Faith's PD)")+
      
      
      ggtitle("Phylogenetic diversity:\nMetabolism")


pd_plot_dgcat <- ggplot(pd_fulldf, aes(x=-dgcat,y=pd)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=pd_fitsdf, aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate),
                  alpha=1.0, color="black", size=2) +
      
      geom_point(size=5, aes(color=Metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +
      #scale_color_manual(values = c_src_colours, labels = c_src_labels, name="Anabolic C source" ) +
      #scale_color_manual(labels = c(C="Complex", S = "Simple"), name="Anabolic C source" ) +
      #xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      xlab(expression(-Delta*G[cat] (kJmol^{-1}))) +
      ylab("Phylogenetic diversity (Faith's PD)")+
      
      ggtitle("Phylogenetic diversity:\nMetabolism") 

ggsave(paste('figures/figure_S1.png', sep=""),pd_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave(paste('figures/figure_1.png', sep=""),pd_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 

#Since there aren't any results for SRB_H2 in the studies examined we remove them from the metacommunity analysis
metabolic_energies <- metabolic_energies[-c(which(metabolic_energies$Metabolism == "SRB_H2")),] 

#Read in the results from the NMGS runs 
metadf = read.table('./Results//metacommunity_diversity_20200817.txt', sep="\t", stringsAsFactors = F)
#metabolic_energies <- metabolic_energies[-c(which(metabolic_energies$Metabolism == "SRB_H2"))]
metadf <- metadf[,c("dirname", as.character(metabolic_energies$Metabolism))]

md <- melt(metadf, value.name = "theta", variable.name = "metabolism", measure.vars = metabolic_energies$Metabolism, )

theta_df <- merge(md, metabolic_energies, by.x = "metabolism", by.y="Metabolism")
#Drop missing values
theta_df <- theta_df[(!is.na(theta_df$theta)) & (theta_df$theta > 0),]

theta_df

theta_fits <- theta_df %>% group_by(dirname) %>% do(data.frame(
    yield=summary(lm(log(theta) ~ log(yield), data=.))$coefficients,
    dgcat=summary(lm(log(theta) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

theta_fits$param = c("intercept","slope")

theta_fitsdf<- melt(theta_fits, id.vars = c("dirname","param")) %>% dcast(dirname ~ param +variable)

theta_fitsdf

#Set up the plot axes - for log scales
xbreaks_simple = c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100),seq(1000,10000,1000))
xlabs_simple <- c("0.001",rep("",9),"0.01",rep("",9),"0.1",rep("",9),"1",rep("",9), "10", rep("",9), "100", "","","","500","","","","", "1000", "", rep("",9))

ybreaks_simple = c(seq(0.1,1,0.1),seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))
ylabs_simple <- c("0.1",rep("",9),"1",rep("",9),"10","20","30","40","50","60","70","80","90","100","",rep("",9),"1000",rep("",9))

theta_df %>% group_by(dirname) %>% filter((dirname=="TARA_ERP001736") & C.Source.type=="S") -> simple_df 
cat("TARA simple carbon fit - yield")
summary(lm(log(theta) ~ log(yield), data=simple_df))
cat("TARA simple carbon fit - DGcat")
summary(lm(log(theta) ~ log(-dgcat), data=simple_df))

theta_df %>% group_by(dirname) %>% filter((dirname=="TARA_ERP001736") & (C.Source.type=="S") &(metabolism != "Methylotroph") ) -> simple_df 
cat("TARA simple carbon fit - yield - dropping Methylotrophs")
summary(lm(log(theta) ~ log(yield), data=simple_df))
cat("TARA simple carbon fit - DGcat - dropping Methylotrophs")
summary(lm(log(theta) ~ log(-dgcat), data=simple_df))

theta_df

theta_fits <- theta_df %>% group_by(dirname) %>% filter((dirname=="TARA_ERP001736") & C.Source.type=="S") %>% do(data.frame(
    yield=summary(lm(log(theta) ~ log(yield), data=.))$coefficients,
    dgcat=summary(lm(log(theta) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

theta_fits$param = c("intercept","slope")
rownames(theta_fitsdf) <- theta_fitsdf$dirname



theta_fitsdf<- melt(theta_fits, id.vars = c("dirname","param")) %>% dcast(dirname ~ param +variable)

theta_plot_yield <- ggplot(theta_df[(theta_df$dirname=="TARA_ERP001736") & (theta_df$C.Source.type=="S"),], aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf[theta_fitsdf$dirname=="TARA_ERP001736",], aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate), alpha=1.0, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks_simple, labels=xlabs_simple) + 
      scale_y_continuous(trans = "log", breaks=ybreaks_simple, labels=ylabs_simple, limits = c(5,30)) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nTARA dataset- simple carbon substrates") 

theta_plot_dgcat <- ggplot(theta_df[(theta_df$dirname=="TARA_ERP001736") & (theta_df$C.Source.type=="S"),], aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf[theta_fitsdf$dirname=="TARA_ERP001736",], aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate), alpha=1.0, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks_simple, labels=xlabs_simple) + 
      scale_y_continuous(trans = "log", breaks=ybreaks_simple, labels=ylabs_simple, limits = c(5,30)) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression(-Delta*G[cat](kJmol^-1))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nTARA dataset - simple carbon substrates") 

ggsave('figures/figure_3.png',theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave('figures//figure_S3.png',theta_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 

theta_plot_dgcat

theta_plot_yield

## Overall fit
biome_theta_df <- theta_df[(theta_df$dirname=="TARA_ERP001736") & (theta_df$C.Source.type=="S"),]
biom_yield = lm(log(theta) ~ log(yield), biome_theta_df)
biom_dgcat = lm(log(theta) ~ log(-dgcat), biome_theta_df)
    
    


summary(biom_yield)


biome_theta_df = merge(theta_df, study_info_df, by.x="dirname", by.y="study")

biome_theta_fits <- biome_theta_df %>% group_by(biome) %>% do(data.frame(
    yield=summary(lm(log(theta) ~ log(yield), data=.))$coefficients,
    dgcat=summary(lm(log(theta) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

biome_theta_fits$param = c("intercept","slope")
#rownames(biome_theta_fits) <- biome_theta_fits$biome

biome_theta_fitsdf<- melt(biome_theta_fits, id.vars = c("biome","param")) %>% dcast(biome ~ param +variable)

biome_theta_fitsdf[c("biome","slope_yield.Estimate","slope_yield.Std..Error", "slope_yield.Pr...t..")]

mdf_casted <- dcast(biome_theta_df[,c("metabolism", "theta", "biome")], metabolism ~ biome, value.var = "theta", fun.aggregate = mean)
mdf_casted2 <-dcast(biome_theta_df[,c("metabolism", "theta", "dirname")] , metabolism ~ dirname, value.var = "theta")


library(Skillings.Mack)

sk_res <- Ski.Mack(as.matrix(1/mdf_casted[,-1]),simulate.p.value = T, B=100 )

sk_res

biome_theta_fitsdf

theta_plot_yield <- ggplot(biome_theta_df, aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=biome_theta_fitsdf, aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate, color=biome),
                  alpha=1.0, size=2) +
      
      #geom_point(size=5, aes(color=biome), alpha=1.0) +
      geom_point(size=2, color="black", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = biome_color_list, name="Source biome" ) +

      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies by biome") 

theta_plot_dgcat <- ggplot(biome_theta_df, aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=biome_theta_fitsdf, aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate, color=biome),
                  alpha=1.0, size=2) +
      
      #geom_point(size=5, aes(color=biome), alpha=1.0) +
      geom_point(size=2, color="black", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = biome_color_list, name="Source biome" ) +

      xlab(expression(-Delta*G[cat](kJmol^-1)))  +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies by biome") 


ggsave('figures/figure_2.png',theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave('figures/figure_S2.png',theta_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
theta_plot_yield
theta_plot_dgcat


anova(lm(log(theta) ~ log(yield)*biome, biome_theta_df))


biome_theta_fitsdf[,c(1,2,5,10,13)][c(1,4,5,2,3)]

summary(lm(log(theta) ~ log(yield):biome + biome + yield , biome_theta_df))

summary(lm(log(theta) ~ log(yield) , biome_theta_df))

summary(lm(log(theta) ~ log(-dgcat) , biome_theta_df))

summary(lm(log(theta) ~ log(yield):biome , biome_theta_df))
