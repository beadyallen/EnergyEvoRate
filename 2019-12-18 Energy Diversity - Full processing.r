
# Set working directory
setwd('/opt/scratch/EnergyDiversity/')

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
biomlist <- list.files('biom_files/', '*.json', full.names = T)
prj_names <- gsub('\\.biom\\.json','\\1',basename(biomlist))
nmgsdirs <- paste('./NMGS_processing/',prj_names, sep="")

prj_names

#Only required if performaing NMGS theta calculations
for (p in nmgsdirs)
#Create the base NMGS directory
    dir.create(p, recursive = T)


classify_data <- read.csv('./Taxa_Classifications_treeprocessed_20191218.csv', header = T, comment.char = '#', stringsAsFactors = F )
metabolic_energies <- read.csv('./Energies_with_RET.csv')

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

PSobjs_rare <- sapply(names(PSobjs_full), function(n){
    ps <- PSobjs_full[[n]]
    h = hist(sample_sums(ps), breaks=40, plot=F)
    top_80 = which(cumsum(h$counts) / sum(h$counts) > 0.2)[[1]]
    rarefy_even_depth(ps, sample.size = h$mids[[top_80]])
    })


tax_table(PSobjs_rare$`biom_files//AmericanGut_ERP012803.biom.json`)

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

#Named list for all results
all_filtered_rare = list()

for (i in seq(length(PSobjs_rare))){
    curbiom = PSobjs_rare[[i]]
    filtered_list <- lapply(Taxa_Search_Expressions, function(x){
                                return(tryCatch(do.call(subset_taxa, list(curbiom,x)),error=function(e) NULL))
                            })
    filtered_list = filtered_list[!sapply(filtered_list, is.null)]
    all_filtered_rare[[prj_names[[i]]]] <- filtered_list
                        
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

apply(as.data.frame(allprojects), 1,mean)

apply(as.data.frame(genusprojects), 1,mean)

for (i in seq(length(all_filtered_list))){
    curbiom <- all_filtered_list[[i]]
    for (x in names(curbiom)){
        outfname = paste(nmgsdirs[[i]],'/',x,'.csv', sep = '')
        write.csv(curbiom[[x]]@otu_table, outfname)
    }
}

read.table('./NMGS_processing/metacommunity_diversity_20180703.txt', sep="\t")

#Set up the plot axes - for log scales
xbreaks = c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100),seq(1000,10000,1000))
xlabs <- c("0.001",rep("",9),"0.01",rep("",9),"0.1",rep("",9),"1",rep("",9), "10", rep("",9), "100", rep("",9), "1000", rep("",9))

ybreaks = c(seq(0.1,1,0.1),seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))
ylabs <- c("0.1",rep("",9),"1",rep("",9),"10",rep("",9),"100",rep("",9),"1000",rep("",9))

#Metabolic colours
metab_color_list <-  c("wheat3", "blue","red", "chartreuse3", "darkgoldenrod2","orange","darkred", "purple", "darkgreen", "yellow4", "steelblue3","firebrick2", "indianred", "darkolivegreen4", "honeydew4")
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




alpha_diversity_outputs <- './output_images/alpha'


#Marine_ERP013553 doesn't have enough diversity per sample (coverage is too low)
all_filtered_list[["Marine_ERP013553"]] <- NULL

all_filtered_rare$AmericanGut_ERP012803

#all_alpha_met_df_list <- list()

alpha_df_list <- lapply(names(all_filtered_list), function(curname){
    curPrj <- all_filtered_list[[curname]]
    
    prj_diversity_df <- do.call(rbind,
            lapply(names(curPrj), function(b){
                df <- estimate_richness(curPrj[[b]], measures = c("Observed","InvSimpson","Shannon", "Chao1"))
                df$metabolism <- b
                df$sample <- rownames(df)
                df$study <- curname
                df
            }))
    prj_diversity_df
})

#all_alpha_met_df_list <- list()

alpha_df_list <- lapply(names(all_filtered_rare), function(curname){
    curPrj <- all_filtered_rare[[curname]]
    
    prj_diversity_df <- do.call(rbind,
            lapply(names(curPrj), function(b){
                df <- estimate_richness(curPrj[[b]], measures = c("Observed","InvSimpson","Shannon", "Chao1"))
                df$metabolism <- b
                df$sample <- rownames(df)
                df$study <- curname
                df
            }))
    prj_diversity_df
})

#combine all diversities into a single data.frame
alpha_diversity = do.call(rbind,alpha_df_list)


#fixup columns 
#Drop "inf" values - not all samples from all projects contain examples from all metabolisms
alpha_diversity <- alpha_diversity[is.finite(alpha_diversity$InvSimpson),]

#and merge in the metabolic values:
alpha_diversity = merge(alpha_diversity, metabolic_energies, by.x="metabolism", by.y="Metabolism")
      

alpha_diversity %>% head

colnames(alpha_diversity)[colnames(alpha_diversity) == "sample"] <- "samp_name"

alpha_diversity  %>% select(metabolism, InvSimpson, samp_name, study) %>% 
    group_by(samp_name) %>% arrange(desc(InvSimpson)) %>% #summarise(metab=list(metabolism), )
    mutate(ord = InvSimpson / max(InvSimpson)) %>% group_by(metabolism) -> tmpdf

tmpdf %>% group_by(metabolism) %>% summarise(m = mean(ord)) %>% arrange(desc(m))

#Calulate per-sample linear best fits
alpha_fits_yield <- alpha_diversity %>% group_by(study, sample) %>% do(
    as.data.frame(summary(lm(log(InvSimpson) ~ log(yield), data=.))$coefficients, stringsAsFactors = F))


alpha_fits_dgcat <- alpha_diversity %>% group_by(study, sample) %>% do(
    as.data.frame(summary(lm(log(InvSimpson) ~ log(-dgcat), data=.))$coefficients, stringsAsFactors = F))

#Finally, group all samples together and generate the "full" fit across tall samples in a project
alpha_fullfits_yield <- alpha_diversity %>% group_by(study) %>% do(as.data.frame(summary(lm(log(InvSimpson) ~ log(yield), data=.))$coefficients, stringsAsFactors = F))
#Finally, group all samples together and generate the "full" fit across tall samples in a project
alpha_fullfits_dgcat <- alpha_diversity %>% group_by(study) %>% do(as.data.frame(summary(lm(log(InvSimpson) ~ log(-dgcat), data=.))$coefficients, stringsAsFactors = F))

#Relabel the parameters (intercept and slope) that were lost "do" aglomeration
alpha_fits_yield$param = c("intercept","slope")
alpha_fullfits_yield$param = c("intercept","slope")
alpha_fits_dgcat$param = c("intercept","slope")
alpha_fullfits_dgcat$param = c("intercept","slope")

#Reshape the dataframes
alpha_fitsdf_yield <- melt(alpha_fits_yield, id.vars = c("study","sample","param")) %>% dcast(study + sample ~ param+variable, value.var = "value")
alpha_fullfitsdf_yield <- melt(alpha_fullfits_yield, id.vars = c("study","param")) %>% dcast(study ~ param+variable, value.var = "value")

alpha_fitsdf_dgcat <- melt(alpha_fits_dgcat, id.vars = c("study","sample","param")) %>% dcast(study + sample ~ param+variable, value.var = "value")
alpha_fullfitsdf_dgcat <- melt(alpha_fullfits_dgcat, id.vars = c("study","param")) %>% dcast(study ~ param+variable, value.var = "value")

rownames(alpha_fullfitsdf_yield) <- alpha_fullfitsdf_yield$study
rownames(alpha_fullfitsdf_dgcat) <- alpha_fullfitsdf_dgcat$study

alpha_fullfits_yield

alpha_fullfits_yield[alpha_fullfits_yield$param == "slope",]

for (curstudy in unique(alpha_fitsdf_dgcat$study)){

    plot_df <- alpha_diversity[alpha_diversity$study == curstudy,]
    
    fits_df_yield <- alpha_fitsdf_yield[alpha_fitsdf_yield$study == curstudy,]
    fits_df_dgcat <- alpha_fitsdf_dgcat[alpha_fitsdf_dgcat$study == curstudy,]



    p_yield <- ggplot(plot_df, aes(x=yield, y=InvSimpson)) +
        theme_bw(base_size = 18, base_family = "arial") +

        #geom_abline(data=fits_df_yield[fits_df_yield$`slope_Pr(>|t|)` <= 0.05,], 
        #              aes(slope=slope_Estimate, intercept=intercept_Estimate), alpha=0.05, color="black") +
        #geom_abline(data=fitdf[fitdf$`slope_Pr(>|t|)` > 0.05,], 
        #              aes(slope=slope_Estimate, intercept=intercept_Estimate), alpha=0.05, color="black") +


        geom_abline(data=fits_df_yield, aes(slope=slope_Estimate, intercept=intercept_Estimate), alpha=0.05, color="black") +

        geom_point(size=5, aes(color=metabolism), alpha=1.0) +
        geom_point(size=2, color="LightBlue", alpha=0.8) + 

        scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
        scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

        scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

        xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
        ylab(expression(alpha ~ "Diversity (inverse Simpson)"))  +
        geom_abline(slope=alpha_fullfits_yield[curstudy, "slope_Estimate"], intercept = alpha_fullfits_yield[curstudy, "intercept_Estimate"],
                 alpha = 1.0, color="darkgreen", size=2) +
        ggtitle(study_name_list[curstudy]) 


    p_catab <- ggplot(plot_df, aes(x=-dgcat, y=InvSimpson)) +
        theme_bw(base_size = 18, base_family = "arial") +

        #geom_abline(data=fits_df_yield[fits_df_yield$`slope_Pr(>|t|)` <= 0.05,], 
        #              aes(slope=slope_Estimate, intercept=intercept_Estimate), alpha=0.05, color="black") +
        #geom_abline(data=fitdf[fitdf$`slope_Pr(>|t|)` > 0.05,], 
        #              aes(slope=slope_Estimate, intercept=intercept_Estimate), alpha=0.05, color="black") +


        geom_abline(data=fits_df_dgcat, aes(slope=slope_Estimate, intercept=intercept_Estimate), alpha=0.05, color="black") +

        geom_point(size=5, aes(color=metabolism), alpha=1.0) +
        geom_point(size=2, color="LightBlue", alpha=0.8) + 

        scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
        scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

        scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +
        xlab(expression(-Delta*G[cat](kJmol^-1))) +

        ylab(expression(alpha ~ "Diversity (inverse Simpson)"))  +
        geom_abline(slope=alpha_fullfits_dgcat[curstudy, "slope_Estimate"], intercept = alpha_fullfits_dgcat[curstudy, "intercept_Estimate"],
                 alpha = 1.0, color="darkgreen", size=2) +
        ggtitle(study_name_list[curstudy]) 


    ggsave(paste('output_images2/alpha/yield-',curstudy,'.png', sep=""),p_yield,dpi = 600, 
               units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 

    ggsave(paste('output_images2/alpha/catab-',curstudy,'.png', sep=""),p_catab,dpi = 600, 
               units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 

}

alpha_fullfits_yield

fits_df_dgcat

#Read in the results from the NMGS runs 
metadf = read.table('./NMGS_processing/metacommunity_diversity_20180703.txt', sep="\t", stringsAsFactors = F)

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

theta_fits

#rownames(theta_fitsdf_yield) <- theta_fitsdf_yield$dirname

#theta_fitsdf<- melt(theta_fits, id.vars = c("dirname","param")) %>% dcast(dirname ~ param +variable)

theta_plot_yield <- ggplot(theta_df, aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf, aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate), alpha=0.5, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies") 

theta_plot_dgcat <- ggplot(theta_df, aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf, aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate), alpha=0.5, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression(-Delta*G[cat](kJmol^-1))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies") 

ggsave('output_images/theta/yield-allstudies.png',theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave('output_images/theta/dgcat-allstudies.png',theta_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 



theta_plot_yield <- ggplot(theta_df[theta_df$dirname == "TARA_ERP001736",], aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf[theta_fitsdf$dirname=="TARA_ERP001736",], aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate), alpha=1.0, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies") 

theta_plot_dgcat <- ggplot(theta_df[theta_df$dirname == "TARA_ERP001736",], aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf[theta_fitsdf$dirname=="TARA_ERP001736",], aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate), alpha=1.0, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression(-Delta*G[cat](kJmol^-1))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies") 

ggsave('output_images/theta/TARA-yield.png',theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave('output_images/theta/TARA-dgcat.png',theta_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 



theta_plot_yield <- ggplot(theta_df, aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf[theta_fitsdf$dirname=="TARA_ERP001736",], aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate), alpha=1.0, color="black") +
      
    
      geom_point(size=5, aes(color=metabolism), alpha=0.1) +
      #geom_point(data=theta_df[theta_df$dirname == "TARA_ERP001736",], aes(x=yield,y=theta), color="black", size=5) +
    
      geom_point(size=2, color="LightBlue", alpha=0.1) + 

      #geom_point(data=theta_df[theta_df$dirname == "TARA_ERP001736",], aes(x=yield,y=theta), color="black", size=5) +
      geom_point(data=theta_df[theta_df$dirname == "TARA_ERP001736",], aes(x=yield,y=theta, color=metabolism), size=5) +
      geom_point(data=theta_df[theta_df$dirname == "TARA_ERP001736",], aes(x=yield,y=theta), color="white", size=3) +
      geom_point(data=theta_df[theta_df$dirname == "TARA_ERP001736",], aes(x=yield,y=theta, color=metabolism), size=1) +
      
      #geom_point(data=theta_df[theta_df$dirname == "TARA_ERP001736",], aes(x=yield,y=theta), color="red", size=2) +

      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies") 

theta_plot_dgcat <- ggplot(theta_df[theta_df$dirname == "TARA_ERP001736",], aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf[theta_fitsdf$dirname=="TARA_ERP001736",], aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate), alpha=1.0, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression(-Delta*G[cat](kJmol^-1))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies") 

ggsave('output_images/theta/TARA-full_yield_highlighted.png',theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave('output_images/theta/TARA-full_dgcat.png',theta_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 



#Set up the plot axes - for log scales
xbreaks = c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100),seq(1000,10000,1000))
xlabs <- c("0.001",rep("",9),"0.01",rep("",9),"0.1",rep("",9),"1",rep("",9), "10", rep("",9), "100", "","","","500","","","","", "1000", "", rep("",9))

ybreaks = c(seq(0.1,1,0.1),seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))
ylabs <- c("0.1",rep("",9),"1",rep("",9),"10","20","30","40","50","60","70","80","90","100","",rep("",9),"1000",rep("",9))

theta_fits <- theta_df %>% group_by(dirname) %>% filter(C.Source.type=="S") %>% do(data.frame(
    yield=summary(lm(log(theta) ~ log(yield), data=.))$coefficients,
    dgcat=summary(lm(log(theta) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

theta_fits$param = c("intercept","slope")
rownames(theta_fitsdf_yield) <- theta_fitsdf_yield$dirname




theta_fitsdf<- melt(theta_fits, id.vars = c("dirname","param")) %>% dcast(dirname ~ param +variable)

theta_plot_yield <- ggplot(theta_df[(theta_df$dirname=="TARA_ERP001736") & (theta_df$C.Source.type=="S"),], aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf[theta_fitsdf$dirname=="TARA_ERP001736",], aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate), alpha=1.0, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs, limits = c(5,30)) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nTARA dataset- simple carbon substrates") 

theta_plot_dgcat <- ggplot(theta_df[(theta_df$dirname=="TARA_ERP001736") & (theta_df$C.Source.type=="S"),], aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf[theta_fitsdf$dirname=="TARA_ERP001736",], aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate), alpha=1.0, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs, limits = c(5,30)) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression(-Delta*G[cat](kJmol^-1))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nTARA dataset - simple carbon substrates") 

ggsave('output_images/theta/simple_carbon_TARA-yield.png',theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave('output_images/theta/simple_carbon_TARA-dgcat.png',theta_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 

theta_plot_dgcat

theta_plot_yield

theta_plot_yield + geom_text(x=1,y=10,label="blah")

?geom_text

theta_fitsdf[theta_fitsdf$dirname=="TARA_ERP001736",]

biome_theta_df = merge(theta_df, study_info_df, by.x="dirname", by.y="study")

biome_theta_fits <- biome_theta_df %>% group_by(biome) %>% do(data.frame(
    yield=summary(lm(log(theta) ~ log(yield), data=.))$coefficients,
    dgcat=summary(lm(log(theta) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

biome_theta_fits$param = c("intercept","slope")
#rownames(biome_theta_fits) <- biome_theta_fits$biome

biome_theta_fitsdf<- melt(biome_theta_fits, id.vars = c("biome","param")) %>% dcast(biome ~ param +variable)

biome_theta_fitsdf[c("biome","slope_yield.Estimate","slope_yield.Std..Error", "slope_dgcat.Pr...t..")]

colnames(biome_theta_fitsdf)

mdf_casted <- dcast(biome_theta_df[,c("metabolism", "theta", "biome")], metabolism ~ biome, value.var = "theta", fun.aggregate = mean)

mdf_casted2 <-dcast(biome_theta_df[,c("metabolism", "theta", "dirname")] , metabolism ~ dirname, value.var = "theta")

mdf_casted2

install.packages("Skillings.Mack")

library(Skillings.Mack)

t(as.matrix(mdf_casted2[,-1]))

sk_res <- Ski.Mack(as.matrix(1/mdf_casted[,-1]),simulate.p.value = T, B=100 )

mdf_casted

t(as.character(mdf_casted$metabolism[order(sk_res$rankdata[,3])]))

sk_res$rankdata[order(sk_res$rankdata[,3]),]

sk_res

as.matrix(mdf_casted[,-1])

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
theta_plot_yield

ggsave('./output_images//theta/yield_allbiomes.png', dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 

anova(lm(log(theta) ~ log(yield)*biome, biome_theta_df))

biome_theta_fitsdf[,c(1,2,5,10,13)][c(1,4,5,2,3)]

library(nlme)

summary(lm(log(theta) ~ log(yield):biome + biome + yield , biome_theta_df))

biome_theta_df

theta_plot_dgcat <- ggplot(biome_theta_df, aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=biome_theta_fitsdf, aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate),
                  alpha=0.5, color="black") +
      
      geom_point(size=5, aes(color=biome), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = biome_color_list, name="Source biome" ) +
      xlab(expression(-Delta*G[cat](kJmol^-1))) +
      
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies by biome") 

#ggsave('output_images/theta/yield-allbiomes.png',theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
#ggsave('output_images/theta/dgcat-allbiomes.png',theta_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 

#theta_plot_yield
#theta_plot_dgcat

theta_fitsdf

for (curbiome in unique(biome_theta_df$biome)){
    
    
    curthetadf = biome_theta_df[biome_theta_df$biome == curbiome,]
    
    sub_theta_df = curthetadf
    
    sub_theta_fits <- sub_theta_df %>% group_by(biome) %>% do(data.frame(
        yield=summary(lm(log(theta) ~ log(yield), data=.))$coefficients,
        dgcat=summary(lm(log(theta) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

    sub_theta_fits$param = c("intercept","slope")

    sub_theta_fitsdf<- melt(sub_theta_fits, id.vars = c("biome","param")) %>% dcast(biome ~ param +variable)

    curthetadf = biome_theta_df[biome_theta_df$biome == curbiome,]
    curfitdf = theta_fitsdf[theta_fitsdf$dirname %in% sample_info_df[sample_info_df$biome == curbiome,"study"], ]
    
    theta_plot_yield <- ggplot(curthetadf, aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +

 
      geom_abline(data=curfitdf, aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate),
                  alpha=0.5, color="black") +
    
      geom_abline(data=sub_theta_fitsdf, aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate),
                  alpha=0.9, color="red", size=1) +
      
      geom_point(size=5, aes(color=biome), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs, limits = c(0.1,1000)) +  

      scale_color_manual(values = biome_color_list, name="Source biome" ) +

      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
    
      ggtitle(paste("Metacommunity diversity:\n", curbiome, sep="\t")) +
      annotate("text", x=0.01, y=1000, label = "blah")
    
          ggsave(paste('output_images/theta/yield-',curbiome,'.png', sep=""),theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE)
    
 }

for (curbiome in unique(biome_theta_df$biome)){
    
    
    curthetadf = biome_theta_df[biome_theta_df$biome == curbiome,]
    sub_theta_df = curthetadf
    
    sub_theta_fits <- sub_theta_df %>% group_by(biome) %>% do(data.frame(
        yield=summary(lm(log(theta) ~ log(yield), data=.))$coefficients,
        dgcat=summary(lm(log(theta) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

    sub_theta_fits$param = c("intercept","slope")

    sub_theta_fitsdf<- melt(sub_theta_fits, id.vars = c("biome","param")) %>% dcast(biome ~ param +variable)

    curthetadf = biome_theta_df[biome_theta_df$biome == curbiome,]
    curfitdf = theta_fitsdf[theta_fitsdf$dirname %in% sample_info_df[sample_info_df$biome == curbiome,"study"], ]
    
    theta_plot_yield <- ggplot(curthetadf, aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +

 
      geom_abline(data=curfitdf, aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate),
                  alpha=0.5, color="black") +
    
      geom_abline(data=sub_theta_fitsdf, aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate),
                  alpha=0.9, color="red", size=1) +
      
      geom_point(size=5, aes(color=biome), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs, limits = c(0.1,2000)) +  

      scale_color_manual(values = biome_color_list, name="Source biome" ) +

      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
    
      ggtitle(paste("Metacommunity diversity:\n", curbiome, sep="\t")) 
    
    
    
   theta_plot_dgcat <- ggplot(curthetadf, aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=curfitdf, aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate),
                  alpha=0.5, color="black") +
      geom_abline(data=sub_theta_fitsdf, aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate),
                  alpha=0.9, color="red", size=1) +
      
      geom_point(size=5, aes(color=biome), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = biome_color_list, name="Source biome" ) +

      
      xlab(expression(-Delta*G[cat](kJmol^-1))) +
      
      ylab(expression(theta ~ "(metacommunity diversity)"))+
    
      ggtitle(paste("Metacommunity diversity:\n", curbiome, sep="\t")) 
    
      ggsave(paste('output_images/theta/dgcat-',curbiome,'.png', sep=""),theta_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
      ggsave(paste('output_images/theta/yield-',curbiome,'.png', sep=""),theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
}

biome_theta_df = merge(theta_df, sample_info_df, by.x="dirname", by.y="study")

sub_theta_df = biome_theta_df[biome_theta_df$biome == "Human",]

sub_theta_fits <- sub_theta_df %>% group_by(biome) %>% do(data.frame(
    yield=summary(lm(log(theta) ~ log(yield), data=.))$coefficients,
    dgcat=summary(lm(log(theta) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

sub_theta_fits$param = c("intercept","slope")

sub_theta_fitsdf<- melt(sub_theta_fits, id.vars = c("biome","param")) %>% dcast(biome ~ param +variable)

sub_theta_fitsdf

sub_theta_fits

theta_plot_dgcat <- ggplot(biome_theta_df, aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=biome_theta_fitsdf, aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate),
                  alpha=0.5, color="black") +
      
      geom_point(size=5, aes(color=biome), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = biome_color_list, name="Source biome" ) +
     
      
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies by biome") 

curbiome="WWT"

sample_info_df[sample_info_df$biome == curbiome,"study"]

theta_plot_dgcat <- ggplot(theta_df, aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=theta_fitsdf, aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate), alpha=0.5, color="black") +
      
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +

      xlab(expression(-Delta*G[cat](kJmol^-1))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))+
      
      
      ggtitle("Metacommunity diversity:\nall studies") 

ggsave('output_images/theta/yield-allstudies.png',theta_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave('output_images/theta/dgcat-allstudies.png',theta_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 

dg_xbreaks = c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000) )
dg_xlabs <- c("0.001",rep("",9),"0.01",rep("",9),"0.1",rep("",9),"1",rep("",9), "10", rep("",9),"100", rep("",9),"1000", rep("",9))

(theta_human_plot <- ggplot(theta_df[theta_df$theta > 0,][theta_df$dirname %in% c("AmericanGut_ERP012803","MalawianChildren_ERP005437"),], aes(x=-dgcat,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
      geom_abline(data=theta_fitdf[c("AmericanGut_ERP012803","MalawianChildren_ERP005437"),], 
                  aes(slope=slope_estimate, intercept=intercept_estimate), size=2,alpha=0.5, color="black") + 
      geom_point(size=5, aes(color=metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_log10(breaks=dg_xbreaks, labels=dg_xlabs) + 
      scale_y_log10(breaks=ybreaks, labels=ylabs) +  
      scale_color_manual(values = color_list) +
      xlab(expression(-Delta*G[cat] (kJmol^{-1}))) +
      
      ylab(expression(theta ~ "(metacommunity diversity)"))  +
      
      ggtitle("Metacommunity diversity\nHuman gut only") 
     )

ggsave("theta_diversity_human_vs_dg.png", plot=theta_energy_plot, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)

theta_fitdf[c("AmericanGut_ERP012803","MalawianChildren_ERP005437"),]

#Load the tree
gg_tree = read_tree_greengenes('gg_13_5_otus_99_annotated.tree')
knowntips = gg_tree$tip.label

#load the database entries
gg_dat = read.table('gg_13_5_taxonomy.txt', sep="\t", header=F, stringsAsFactors = F)
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

dim(cleandf)

ul = unique(unlist(tax_groups))
m = match(gg_tree$tip.label, ul, )
prunelist = gg_tree$tip.label[is.na(m)]

nvec=names(tax_groups)
names(nvec)=nvec
all_prune_taxa = lapply(nvec, function(x) unique(unlist(tax_groups[-which(nvec == x)])))

all_prune_taxa

working_tree = drop.tip(gg_tree, prunelist)

all_trees = lapply(all_prune_taxa, function(x) drop.tip(working_tree,x))

pd_list = lapply(all_trees, function(x) sum(x$edge.length))
pd_df = as.data.frame(t(data.frame(pd_list)))
pd_df[,2] = rownames(pd_df)
colnames(pd_df) <- c("pd","metab")

fulldf = merge(metabolic_energies, pd_df, by.x="Metabolism", by.y="metab")
pd_fulldf <- fulldf[order(fulldf$yield, decreasing = T),]

pd_fulldf 

pd_fitsdf

pd_fits <- pd_fulldf %>% filter(pd > 0.0) %>% do(data.frame(
        yield=summary(lm(log(pd) ~ log(yield), data=.))$coefficients,
        dgcat=summary(lm(log(pd) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

pd_fits$param = c("intercept","slope")

pd_fitsdf<- melt(pd_fits, id.vars = c("param")) %>% dcast(. ~ param +variable)



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

ggsave(paste('output_images/phylogeny/PD-dgcat-metabolisms-.png', sep=""),pd_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave(paste('output_images/phylogeny/PD-yield-metabolisms-.png', sep=""),pd_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 


pd_fits <- pd_fulldf %>% filter(pd > 0.0) %>% group_by(C.Source.type) %>% do(data.frame(
        yield=summary(lm(log(pd) ~ log(yield), data=.))$coefficients,
        dgcat=summary(lm(log(pd) ~ log(-dgcat), data=.))$coefficients,stringsAsFactors = F))

pd_fits$param = c("intercept","slope")

pd_fitsdf<- melt(pd_fits, id.vars = c("C.Source.type","param")) %>% dcast(C.Source.type ~ param +variable)



pd_plot_yield <- ggplot(pd_fulldf, aes(x=yield,y=pd)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=pd_fitsdf, aes(slope=slope_yield.Estimate, intercept=intercept_yield.Estimate), alpha=0.5, color="black") +
      
      geom_point(size=5, aes(color=C.Source.type), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      #scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +
      #scale_color_manual(values = metab_color_list, labels = metab_label_list, name="Metabolic classification" ) +
      scale_color_manual(values = c_src_colours, labels = c_src_labels, name="Anabolic C source" ) +
      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab("Phylogenetic diversity (Faith's PD)")+
      
      
      ggtitle("Phylogenetic diversity:\nCarbon source") 


pd_plot_dgcat <- ggplot(pd_fulldf, aes(x=-dgcat,y=pd)) + theme_bw(base_size = 18, base_family = "arial") +
 
      geom_abline(data=pd_fitsdf, aes(slope=slope_dgcat.Estimate, intercept=intercept_dgcat.Estimate), alpha=0.5, color="black") +
      
      geom_point(size=5, aes(color=C.Source.type), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_continuous(trans = "log", breaks=xbreaks, labels=xlabs) + 
      scale_y_continuous(trans = "log", breaks=ybreaks, labels=ylabs) +  

      scale_color_manual(values = c_src_colours, labels = c_src_labels, name="Anabolic C source" ) +
      #scale_color_manual(labels = c(C="Complex", S = "Simple"), name="Anabolic C source" ) +
      #xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      xlab(expression(-Delta*G[cat] (kJmol^{-1}))) +
      ylab("Phylogenetic diversity (Faith's PD)")+
      
      
      ggtitle("Phylogenetic diversity:\nCarbon source") 


ggsave(paste('output_images/phylogeny/PD-dgcat-c-src-.png', sep=""),pd_plot_dgcat,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 
ggsave(paste('output_images/phylogeny/PD-yield-c-src-.png', sep=""),pd_plot_yield,dpi = 600,units="mm", width = 150*2, height = 100*2, limitsize = FALSE) 







pd_fulldf

!pd_fulldf$Metabolism %in% c("AerobicResp","Ferm","SRB", "ANAMOX", "SRB_S0")

pd_fulldf <- fulldf[order(fulldf$yield, decreasing = T),]

pd_fulldf <- pd_fulldf[,!pd_fulldf$Metabolism %in% c("AerobicResp","Ferm","SRB", "ANAMOX", "SRB_S0")]
pdfit <- summary(lm(log10(pd) ~ log10(yield), pd_fulldf))$coefficients


pd_plot <- ggplot(pd_fulldf, aes(x=yield,y=pd)) + theme_bw(base_size = 18, base_family = "arial") +
      geom_point(size=5, aes(color=Metabolism), alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 
      scale_x_log10(breaks=xbreaks, labels=xlabs) + 
      scale_y_log10(breaks=ybreaks, labels=ylabs) +  
      scale_color_manual(values = color_list) +
      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression("PD (Faith's distance)"))  +
      geom_abline(slope=pdfit[["log10(yield)","Estimate"]], intercept = pdfit[["(Intercept)","Estimate"]],
                 alpha = 1.0, color="darkgreen", size=2) + 
      ggtitle("Phylogenetic distance\n\"Simple\" metabolisms") 
     

pd_fulldf

pd_plot

pdfit

ggsave("pd_diversity_vs_yield_NoSRBFermAerobic.png", plot=pd_plot, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)

summary(pdfit)

ggsave("pd_diversity_vs_yield.png", plot=pd_plot, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)

Taxa_Search_Expressions$ANAMOX

PSobjs_full[["biom_files//AmericanGut_ERP012803.biom.json"]]

hist(otu_table(all_filtered_list$AmericanGut_ERP012803$AerobicResp))

d = apply(otu_table(all_filtered_list$AmericanGut_ERP012803$AerobicResp), 2, function(x) hist(x,plot=F, breaks = seq(0,10000,10)))

library(BiodiversityR)

df <- as.data.frame(otu_table(PSobjs_full[["biom_files//AmericanGut_ERP012803.biom.json"]]))
tab = table(rowSums(df[rowSums(df)>1,]))

y=as.vector(tab)
x=as.vector(as.numeric(names(tab)))

df = data.frame("freq"=y, "abund"=x)
dfamer <- df
dfamer$source <- "American gut"

P1 <- ggplot(df,aes(x=abund, y=freq)) + geom_point() +
    scale_y_log10() + 
    xlab(expression("Frequency")) +
    ylab(expression("OTU abundance"))  +
    ggtitle("American Gut Project - overall abundance") +
    theme_bw() + 
    scale_x_continuous(limits = c(0,1000))

df <- as.data.frame(otu_table(PSobjs_full[["biom_files//MalawianChildren_ERP005437.biom.json"]]))

tab = table(rowSums(df[rowSums(df)>1,]))

y=as.vector(tab)
x=as.vector(as.numeric(names(tab)))

df = data.frame("freq"=y, "abund"=x)
df_malawi <- df

df_malawi$source <- "Malwian children"

P2 <- ggplot(df,aes(x=abund, y=freq)) + geom_point() +
    scale_y_log10() + 
    xlab(expression("Frequency")) +
    ylab(expression("OTU abundance"))  +
    ggtitle("Malawian children - overall abundance") +
    theme_bw() + 
    scale_x_continuous(limits = c(0,1000))

dfamer

str(rbind(dfamer,df_malawi))

fulldf <- rbind(dfamer,df_malawi)

tail(fulldf)

P1

P2


P3 <- ggplot(rbind(df_malawi,dfamer),aes(x=abund, y=freq, group=source)) + geom_point(aes(color=source, alpha=0.6)) +
    scale_y_log10() + 
    xlab(expression("Frequency")) +
    ylab(expression("OTU abundance"))  +
    ggtitle("Combined - overall abundance") +
    theme_bw() + 
    scale_x_continuous(limits = c(0,1000))

ggsave("malawianchildren_freqabund.png", plot=P2, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)
ggsave("americangut_freqabund.png", plot=P1, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)

ggsave("combined_human_freqabund.png", plot=P3, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)

hist(df[,1], plot=F)

df[1,]

df <- as.data.frame(otu_table(PSobjs_full[["biom_files//AmericanGut_ERP012803.biom.json"]]))

tab = df[rowSums(df)>1,]

obs_freq = apply(tab,1,function(x) sum(x>0) / length(x) * 100) 

mean_abund = apply(apply(tab,2,function(x) x/sum(x)),1,mean)

pidf <- data.frame("of" = obs_freq, "ab"= mean_abund)

(p1<-ggplot(pidf, aes(x=ab, y=of)) + geom_point() + scale_x_continuous(limits = c(0,0.025)) + 
theme_bw() + 
    scale_x_log10(limits=c(0.00001,1.2), breaks=c(0.0001, 0.001, 0.01,0.1,1), labels=c("0.0001","0.001","0.01","0.1","1")) +
    labs(x="Mean relative abundance", y="Detection frequency") + guides(color=FALSE) +
    ggtitle('American Gut project\nDetection frequency vs mean abundance'))

df <- as.data.frame(otu_table(PSobjs_full[["biom_files//MalawianChildren_ERP005437.biom.json"]]))

tab = df[rowSums(df)>1,]

obs_freq = apply(tab,1,function(x) sum(x>0) / length(x) * 100) 

mean_abund = apply(apply(tab,2,function(x) x/sum(x)),1,mean)

pidf <- data.frame("of" = obs_freq, "ab"= mean_abund)

(p2<-ggplot(pidf, aes(x=ab, y=of)) + geom_point() + scale_x_continuous(limits = c(0,0.025)) + 
theme_bw() + 
    scale_x_log10(limits=c(0.00001,1.2), breaks=c(0.0001, 0.001, 0.01,0.1,1), labels=c("0.0001","0.001","0.01","0.1","1")) +
    labs(x="Mean relative abundance", y="Detection frequency") + guides(color=FALSE) +
    ggtitle('Malawian children\nDetection frequency vs mean abundance'))

ggsave("malawianchildren_freqabund_mean.png", plot=p2, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)
ggsave("americangut_freqabund_mean.png", plot=p1, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)



df <- as.data.frame(otu_table(PSobjs_full[["biom_files//MalawianChildren_ERP005437.biom.json"]]))

tab = df[rowSums(df)>1,]

Bact <- data.frame(t(tab))
#Bacs <- data.frame(t(bac_s@otu_table))

Bact <- Bact[order(rownames(Bact)),]
#Bacs <- Bacs[order(rownames(Bacs)),]

Bacs <- replicate(dim(Bact)[1],apply(tab,1,mean))

Bacs <- Bacs[order(rownames(Bacs)),]

#data.frame("src"=rownames(Bacs), "dst" = rownames(Bact))
allOTUs <- 1:ncol(Bact)


Bactfreq=c()
for (i in 1:ncol(Bact)){
        Bactfreq[i]=length(which(Bact[,i]>0))/nrow(Bact)
        }



Bactpi=colMeans(Bact)



Bacsfreq=c()
for (i in 1:ncol(Bacs)){
        Bacsfreq[i]=length(which(Bacs[,i]>0))/nrow(Bacs)
        }

Bacspi=colMeans(Bacs)

colind=c(1:ncol(Bact))

Bac_neutral_data=matrix(nrow=length(Bactfreq),ncol=6)

for (i in 1:nrow(Bac_neutral_data)){

        Bac_neutral_data[i,1]=Bactfreq[i]
        Bac_neutral_data[i,2]=Bactpi[i]
        Bac_neutral_data[i,3]=Bacsfreq[i]
        Bac_neutral_data[i,4]=Bacspi[i]
        Bac_neutral_data[i,5]=colind[i]
        Bac_neutral_data[i,6]=allOTUs[i]
}


# Sort the OTUs by their abundance in the source community just like the neutral model in the excel spreadsheet for the neutral model 
Bac_neutral_data_sorted=Bac_neutral_data[order(Bac_neutral_data[,4],decreasing=TRUE),]
reactor=length(which(Bac_neutral_data_sorted[,1]>0))
source=length(which(Bac_neutral_data_sorted[,3]>0))

shared=length(which(Bac_neutral_data_sorted[,1]>0 & Bac_neutral_data_sorted[,3]>0))
#quartz()
#draw.pairwise.venn(source,reactor,shared,col=c("red","black"),fill=c("indianred","black"),label.col=c("white"))


# consider only shared OTUs b/w the target and source
onlyshared=which(Bac_neutral_data_sorted[,1]>0)
Bac_neutral_data_final=Bac_neutral_data_sorted[onlyshared,]

# 6.25e-5 is the minimum abundance detected
col=which(Bac_neutral_data_final[,4]>=8.6e-7)
Bac_neutral_data_final_1=matrix(nrow=length(col),ncol=6)
for (i in 1:length(col)){

        Bac_neutral_data_final_1[i,1]= Bac_neutral_data_final[i,1]
        Bac_neutral_data_final_1[i,2]= Bac_neutral_data_final[i,2]
        Bac_neutral_data_final_1[i,3]= Bac_neutral_data_final[i,3]
        Bac_neutral_data_final_1[i,4]= Bac_neutral_data_final[i,4]
        Bac_neutral_data_final_1[i,5]= Bac_neutral_data_final[i,5]
        Bac_neutral_data_final_1[i,6]= Bac_neutral_data_final[i,6]
}

obs <- Bac_neutral_data_final_1

NtmFunc = function(Nt, obs){
    detlim = min(obs[,4])
    
    neutralmatrix= matrix(nrow=nrow(obs),ncol=5)
    
    for (i in 1:nrow(obs)){
            neutralmatrix[i,1]=  Nt*obs[i,4]
        }
    for (i in 1:nrow(obs)){
            neutralmatrix[i,2]=  Nt*(1-obs[i,4])
        }
    for (i in 1:nrow(obs)){
            neutralmatrix[i,3]= pbeta(detlim,neutralmatrix[i,1],neutralmatrix[i,2])
        }
    for (i in 1:nrow(obs)){
            neutralmatrix[i,4]= 1-neutralmatrix[i,3]
        }
    for (i in 1:nrow(obs)){
            neutralmatrix[i,5]= (obs[i,1]-neutralmatrix[i,4])^2
        }
    sumcolumns=colSums(neutralmatrix)
    return (sumcolumns[5])
}

bestNtm<- optim(c(100), function(x) NtmFunc(x, obs), method = "Brent", lower = 0, upper = 100000000000000000000)$par
    
bestNtm

obs_freq = apply(tab,1,function(x) sum(x>0) / length(x) * 100) 

mean_abund = apply(apply(tab,2,function(x) x/sum(x)),1,mean)

pidf <- data.frame("of" = obs_freq, "ab"= mean_abund)

(p2<-ggplot(pidf, aes(x=ab, y=of)) + geom_point() + scale_x_continuous(limits = c(0,0.025)) + 
theme_bw() + 
    scale_x_log10(limits=c(0.00001,1.2), breaks=c(0.0001, 0.001, 0.01,0.1,1), labels=c("0.0001","0.001","0.01","0.1","1")) +
    labs(x="Mean relative abundance", y="Detection frequency") + guides(color=FALSE) +
    ggtitle('Malawian children\nDetection frequency vs mean abundance')

all_filtered_list$AmericanGut_ERP012803$AerobicResp

df <- as.data.frame(otu_table(all_filtered_list$AmericanGut_ERP012803$AerobicResp))

tab = df[rowSums(df)>1,]
tab = df[,colSums(df)>0]

obs_freq = apply(tab,1,function(x) sum(x>0) / length(x) * 100) 

mean_abund = apply(apply(tab,2,function(x) x/sum(x, na.rm = T)),1,mean)

pidf <- data.frame("of" = obs_freq, "ab"= mean_abund)

(p2<-ggplot(pidf, aes(x=ab, y=of)) + geom_point()+ 
theme_bw() + 
    scale_x_log10(limits=c(0.00001,1.2), breaks=c(0.0001, 0.001, 0.01,0.1,1), labels=c("0.0001","0.001","0.01","0.1","1")) +
    labs(x="Mean relative abundance", y="Detection frequency") + guides(color=FALSE) +
    ggtitle('Malawian children\nDetection frequency vs mean abundance'))

library(plyr)
library(dplyr)

aldf_list <- list()
for (curname in names(all_filtered_list)){
    curPrj <- all_filtered_list[[curname]]
    
    alpha_diversity_df <- do.call(rbind,
            lapply(names(curPrj), function(b){
                df <- estimate_richness(curPrj[[b]], measures = c("Observed","InvSimpson","Shannon", "Chao1"))
                df$metabolism <- b
                df$sample <- rownames(df)
                df
            }))
    
    alpha_diversity_df[["project"]] <- curname
    aldf_list[[curname]] <- alpha_diversity_df 
}

aldf <- rbind.fill(aldf_list)

#Drop NA entries
aldf <- aldf[! apply(is.na(aldf), 1, any),]

library(reshape2)


casted_df  <- dcast(aldf, project + sample ~ metabolism, value.var = "InvSimpson" )

subcasted <- casted_df[,c('Acet_Methanogen', 'AerobicResp', 'Ammonia_Oxidizer', 'ANAMOX', 'CH3_Methanogen', 'Fe_oxidizer', 'FeIII_reducer', 'Ferm' , 'H2_Methanogen', 'H2_Oxidizer', 'Methanotroph', 'NOB', 'SOB', 'SRB')]

friedman.test(as.matrix(subcasted), na.action = )

getOption("na.action")

library(Skillings.Mack)


mat <- as.matrix(subcasted[apply(!is.na(subcasted), 1, sum) > 1,])

table(apply(mat >0,1,function(x) sum(x,na.rm = T)))

Ski.Mack(mat)

#Drop single values


Ski.Mack(as.matrix(subcasted[apply(!is.na(subcasted), 1, sum) > 2,]))

str(head(as.matrix(subcasted[apply(!is.na(subcasted), 1, sum) > 2,])))

maty <- matrix(
c(3.2,4.1,3.8,4.2,3.1,3.9,3.4,4.0, 4.3,3.5,4.6,4.8,
3.5,3.6,3.9,4.0, 3.6,4.2,3.7,3.9, 4.5,4.7,3.7, NA,
NA ,4.2,3.4,NA , 4.3,4.6,4.4,4.9, 3.5, NA,NA, 34),
ncol=9,byrow=FALSE)
Ski.Mack(maty)

maty

classify_data <- read.csv('./inputs/Taxa_Classifications.csv', header = T, comment.char = '#', stringsAsFactors = F )
met_paths <- unique(classify_data$Metabolism)
names(met_paths) <- met_paths

#Build up Taxa search expressions for them
Taxa_Search_Expressions <- lapply(met_paths, FUN = function(n){
  X = classify_data[classify_data$Metabolism==n,]
  d <- apply(X, 1, function(X) paste(X[[2]],'=="',X[[1]],'"',sep=''))
  searchstr  = paste(unlist(d),collapse = ' | ')
  parse(text=searchstr)
})

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

t1 = tax_table(all_filtered_list$AmericanGut_ERP012803$AerobicResp)
t2 =  tax_table(all_filtered_list$AmericanGut_ERP012803$Ferm)

unique(t1[t1[,"Genus"] %in% t2[,"Genus"],])

all_met_df <- do.call(rbind, all_alpha_met_df_list)

corlist = lapply(all_met_df$sample, function(cursample) {
#cursample = "ERR844083"

    tmp  = all_met_df[all_met_df$sample==cursample,]
    row.names(tmp) <- tmp$metabolism

    m=matrix(0,nrow = length(tmp$metabolism), ncol=length(tmp$metabolism) ,dimnames=list(tmp$metabolism, tmp$metabolism))

    for (i in tmp$metabolism){
        for (j in tmp$metabolism){
            m[i,j] = tmp$InvSimpson[tmp$metabolism == i] / tmp$InvSimpson[tmp$metabolism == j]
        }

    }

    m[lower.tri(m, diag = T)] <- NA


    #m <- m[upper.tri(m)]
    d <- as.data.frame(m)
    d$metabolism <- row.names(m)
    md <- melt(d, id.vars = c("metabolism"))
    colnames(md) <- c("met1","met2","ratio")
    md <- as.data.frame(na.omit(md))
    if (nrow(md) >0 ){
        md$samplename <- cursample
        md
    }else{NA}
})

at_least_3_entries = corlist[unlist(lapply(corlist, function(x) nrow(x) > 2))]

cor_df <- do.call(rbind, corlist)

at_least_3 <- cor_df %>% group_by(samplename) %>% filter(n() > 2)

theta_df

corlist = lapply(theta_df$dirname, function(cursample) {
#cursample = "ERR844083"

    tmp  = theta_df[theta_df$dirname==cursample,]
    row.names(tmp) <- tmp$metabolism

    m=matrix(0,nrow = length(tmp$metabolism), ncol=length(tmp$metabolism) ,dimnames=list(tmp$metabolism, tmp$metabolism))

    for (i in tmp$metabolism){
        for (j in tmp$metabolism){
            m[i,j] = log10(tmp$theta[tmp$metabolism == i]) / log10(tmp$theta[tmp$metabolism == j])
        }

    }

    m[lower.tri(m, diag = T)] <- NA


    #m <- m[upper.tri(m)]
    d <- as.data.frame(m)
    d$metabolism <- row.names(m)
    md <- melt(d, id.vars = c("metabolism"))
    colnames(md) <- c("met1","met2","ratio")
    md <- as.data.frame(na.omit(md))
    if (nrow(md) >0 ){
        md$samplename <- cursample
        md
    }else{NA}
})

cor_theta_df <- do.call(rbind, corlist)

cor_theta_df$met_comb <- paste(cor_theta_df$met1, cor_theta_df$met2, sep="_")

cor_theta_df$gt = cor_theta_df$ratio > 1

cor_theta_df %>% head

library(lme4)

f0 <- lm(ratio ~ met_comb, data = cor_theta_df)
f1 <- lme(ratio ~ met_comb, random = ~ 1 | samplename, data = cor_theta_df)

summary(f1)

f <- lmer(ratio ~ met_comb +(1|samplename), data = cor_theta_df)

anova(f)

rel_theta_df <- theta_df %>% group_by(dirname) %>% mutate("reltheta" = theta / sum(theta)) %>% mutate("rank" = order(reltheta, decreasing = T))

rel_theta_df %>% group_by(dirname) %>% summarise("numentries" = n())

#Fill in metabolisms that aren't there...
mdf <- rel_theta_df[,c("metabolism", "dirname", "reltheta")]

metablist = unique(mdf$metabolism)


for (d in mdf$dirname){
    for (m in metablist){
        if (! any(mdf$dirname == d & mdf$metabolism == m)){
            mdf <- rbind(as.data.frame(mdf), data.frame("dirname" = c(d), "metabolism" = c(m), "reltheta" = c(0.0)))
        }
    }
}

rel_theta_rank_df <- mdf %>% mutate("rank" = order(reltheta, decreasing = T))

f <- lm(rank ~ metabolism , data=rel_theta_rank_df)
summary(f)

f <- lme(rank ~ metabolism , random = ~1|dirname, data=rel_theta_rank_df)
summary(f)

rel_theta_rank_df <- merge(rel_theta_rank_df, biomdf, by.x="dirname", by.y="samplename")

f <- lme(rank ~ metabolism , random = ~1|biom, data=rel_theta_rank_df)
summary(f)

str(rel_theta_rank_df)

rel_alpha_df <- all_met_df %>% group_by(project) %>% mutate("relalpha" = InvSimpson / sum(InvSimpson)) 

#Fill in metabolisms that aren't there...
mdf <- rel_alpha_df[,c("metabolism", "project", "relalpha")]

metablist = unique(mdf$metabolism)

for (d in mdf$project){
    for (m in metablist){
        if (! any(mdf$project == d & mdf$metabolism == m)){
            mdf <- rbind(as.data.frame(mdf), data.frame("project" = c(d), "metabolism" = c(m), "relalpha" = c(0.0)))
        }
    }
}

rel_alpha_rank_df <- mdf %>% mutate("rank" = order(relalpha, decreasing = T))

f <- lm(rank ~ metabolism , data=rel_alpha_rank_df)
summary(f)

f <- lme(rank ~ metabolism , random = ~1|project, data=rel_alpha_rank_df)
summary(f)

#linear mixed effects model

(lme(log10(theta) ~ log10(yield) , random = ~1|dirname, data = theta_df))

summary(lm(theta~yield, data=theta_df))

f1 <- lm(log10(theta) ~ log10(yield), data=theta_df)
f3 <- lm(log10(theta) ~ log10(yield)/dirname  , data=theta_df)
anova(f1,f3)

#theta_df <- merge(theta_df, biomdf, by.x="dirname", by.y="samplename")
f1 <- lm(log10(theta) ~ log10(yield), data=theta_df)
f3 <- lm(log10(theta) ~ log10(yield)/biom  + biom , data=theta_df)
anova(f1,f3)

summary(f3)

f1 <- lm(log10(InvSimpson) ~ log10(yield) + project, all_met_df)
f3 <- lm(log10(InvSimpson) ~ log10(yield)/project + project, all_met_df)

anova(f1,f3)

summary(f3)

sub_met_df <- all_met_df[all_met_df$project == "Lise_Soils",]
f1 <- lm(log10(InvSimpson) ~ log10(yield), data=sub_met_df)
f3 <- lm(log10(InvSimpson) ~ log10(yield) + sample, data=sub_met_df)
anova(f1,f3)

anova(f3)

unique(all_met_df$project)

colnames(all_met_df)



xy_xx <- theta_df %>% group_by(dirname) %>% summarise(xy = sum(theta*yield), xx = sum(yield*yield))
xy_xx <- merge(xy_xx, theta_fitdf[,c("slope_estimate","samplename")], by.x = "dirname", by.y="samplename")

(cs <- colSums(xy_xx[,c(2,3)]))

(bhat = cs["xy"]/cs["xx"])

bhat

log10(xy_xx$xy) / log10(xy_xx$xx)

xy_xx

theta_fitdf

theta_fitdf$biom <- c("Marine", "Freshwater", "Wastewater", "Wastewater", "Soil", "WasteWater", "Human", "Volcanic", "Human",
                     "Marine", "Soil", "Soil", "Volcanic")

biomdf <- theta_fitdf[,c("samplename","biom")]

biomdf

theta_df %>% head

nlist = unique(theta_df$biom)
theta_fitlist <- lapply(nlist,
                      function(x){
                          tmpdf <- theta_df[theta_df$biom == x,]

                          f <- lm( log10(theta) ~ log10(yield), data = tmpdf)
                          summary(f)$coefficients
                      })
names(theta_fitlist) <- nlist
#Filter out samples with insufficient points

theta_fitlist <- theta_fitlist[unlist(lapply(theta_fitlist, function(x) !any(is.na(x))))]


tmplist <- lapply(names(theta_fitlist), function(x){
        tmpdf <- as.data.frame(theta_fitlist[[x]], stringsAsFactors = F)
        tmpdf$samplename <- x
        tmpdf <- t(unlist(tmpdf))
        tmpdf <- t(tmpdf[,-10])
        colnames(tmpdf) <- c("intercept_estimate", "slope_estimate", "intercept_error", "slope_error",
                         "intercept_tval", "slope_tval","intercept_pval", "slope_pval","samplename")
        tmpdf
    })



theta_fitdf = Reduce(function(...) merge(..., all=T), tmplist)
rownames(theta_fitdf) <- theta_fitdf$samplename

    #fitdf <- as.data.frame(t(as.data.frame(fitlist)))
    #colnames(fitdf) <- c("intercept_estimate", "slope_estimate", "intercept_error", "slope_error",
    #                     "intercept_tval", "slope_tval","intercept_pval", "slope_pval")
    for (n in colnames(theta_fitdf)[1:8]){
        theta_fitdf[,n] <- as.numeric(as.character(theta_fitdf[,n]))
    }

save(theta_df, theta_fitdf, file = "ThetaBiom_data.RData")

#Set up some plotting constants - jsut to look nice...
xbreaks = c(seq(0.001,0.01,0.001),seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1), seq(10,100,10))
xlabs <- c("0.001",rep("",9),"0.01",rep("",9),"0.1",rep("",9),"1",rep("",9), "10", rep("",9))

ybreaks = c(seq(0.1,1,0.1),seq(1,10,1), seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))
ylabs <- c("0.1",rep("",9),"1",rep("",9),"10",rep("",9),"100",rep("",9),"1000",rep("",9))

color_list <-  c("wheat3", "blue","red", "green", "darkgoldenrod2","orange","darkred", "purple", "darkgreen", "yellow4", "steelblue3","firebrick2", "indianred", "darkolivegreen4", "honeydew4", "thistle3")
names(color_list) <- metabolic_energies$Metabolism


theta_fitdf

load("ThetaBiom_data.RData")
theta_plot <- ggplot(theta_df[theta_df$theta > 0,], aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
      #geom_abline(data=theta_fitdf[theta_fitdf$slope_pval <= 0.05,], 
      #            aes(slope=slope_estimate, intercept=intercept_estimate), alpha=0.5, color="blue") +
      geom_abline(data=theta_fitdf, 
                  aes(slope=slope_estimate, intercept=intercept_estimate, color=samplename), alpha=0.8, size=2) + 

      geom_point(size=5, alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 

      scale_x_log10(breaks=xbreaks, labels=xlabs) + 
      scale_y_log10(breaks=ybreaks, labels=ylabs) +  
      #scale_color_manual(values = color_list) +
      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))  +
      
      ggtitle("Metacommunity diversity\nacross biomes") 
     

theta_fitdf

theta_plot

#load("ThetaBiom_data.RData")


theta_df <- merge(md, metabolic_energies, by.x = "metabolism", by.y="Metabolism")
#Drop missing values
theta_df <- theta_df[(!is.na(theta_df$theta)) & (theta_df$theta > 0),]

theta_df <- theta_df[(theta_df$theta > 0) & (!pd_fulldf$Metabolism %in% c("AerobicResp","Ferm","SRB", "ANAMOX", "SRB_S0")) ,]

biomdf

md

theta_df <- merge(md, metabolic_energies, by.x = "metabolism", by.y="Metabolism")
#Drop missing values
theta_df <- theta_df[(!is.na(theta_df$theta)) & (theta_df$theta > 0),]

#theta_df <- theta_df[(theta_df$theta > 0) & (!pd_fulldf$Metabolism %in% c("AerobicResp","Ferm","SRB", "ANAMOX", "SRB_S0")) ,]

theta_df <- merge(theta_df, biomdf, by.x="dirname", by.y="samplename")

#theta_df <- theta_df[!theta_df$metabolism %in% c("AerobicResp","Ferm","SRB", "ANAMOX"),]


theta_df

nlist = unique(theta_df$dirname)
theta_fitlist <- lapply(nlist,
                      function(x){
                          tmpdf <- theta_df[theta_df$dirname == x,]

                          f <- lm( log10(theta) ~ log10(yield), data = tmpdf)
                          summary(f)$coefficients
                      })
names(theta_fitlist) <- nlist
#Filter out samples with insufficient points

theta_fitlist <- theta_fitlist[unlist(lapply(theta_fitlist, function(x) !any(is.na(x))))]


tmplist <- lapply(names(theta_fitlist), function(x){
        tmpdf <- as.data.frame(theta_fitlist[[x]], stringsAsFactors = F)
        tmpdf$samplename <- x
        tmpdf <- t(unlist(tmpdf))
        tmpdf <- t(tmpdf[,-10])
        colnames(tmpdf) <- c("intercept_estimate", "slope_estimate", "intercept_error", "slope_error",
                         "intercept_tval", "slope_tval","intercept_pval", "slope_pval","samplename")
        tmpdf
    })



theta_fitdf = Reduce(function(...) merge(..., all=T), tmplist)
rownames(theta_fitdf) <- theta_fitdf$samplename

    #fitdf <- as.data.frame(t(as.data.frame(fitlist)))
    #colnames(fitdf) <- c("intercept_estimate", "slope_estimate", "intercept_error", "slope_error",
    #                     "intercept_tval", "slope_tval","intercept_pval", "slope_pval")
    for (n in colnames(theta_fitdf)[1:8]){
        theta_fitdf[,n] <- as.numeric(as.character(theta_fitdf[,n]))
    }

theta_fitdf

theta_plot <- ggplot(theta_df, aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
      #geom_abline(data=theta_fitdf[theta_fitdf$slope_pval <= 0.05,], 
      #            aes(slope=slope_estimate, intercept=intercept_estimate), alpha=0.5, color="blue") +
      geom_abline(data=theta_fitdf, 
                  aes(slope=slope_estimate, intercept=intercept_estimate, color=samplename), alpha=0.8, size=2) + 

      geom_point(size=5, alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 

      scale_x_log10(breaks=xbreaks, labels=xlabs) + 
      scale_y_log10(breaks=ybreaks, labels=ylabs) +  
      #scale_color_manual(values = color_list) +
      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))  +
      
      ggtitle("Metacommunity diversity\nacross biomes\n\"Simple\" metabolisms") 
     

theta_plot

ggsave("theta_diversity_vs_yield_by_study_simple_metabolisms.png", plot=theta_plot, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)

theta_fitdf

theta_df <- merge(theta_df, biomdf, by.x="dirname", by.y="samplename")

theta_df <- theta_df[!theta_df$metabolism %in% c("AerobicResp","Ferm","SRB", "ANAMOX"),]

theta_df

nlist = unique(theta_df$biom)
theta_fitlist <- lapply(nlist,
                      function(x){
                          tmpdf <- theta_df[theta_df$biom == x,]

                          f <- lm( log10(theta) ~ log10(yield), data = tmpdf)
                          summary(f)$coefficients
                      })
names(theta_fitlist) <- nlist
#Filter out samples with insufficient points

theta_fitlist <- theta_fitlist[unlist(lapply(theta_fitlist, function(x) !any(is.na(x))))]


tmplist <- lapply(names(theta_fitlist), function(x){
        tmpdf <- as.data.frame(theta_fitlist[[x]], stringsAsFactors = F)
        tmpdf$samplename <- x
        tmpdf <- t(unlist(tmpdf))
        tmpdf <- t(tmpdf[,-10])
        colnames(tmpdf) <- c("intercept_estimate", "slope_estimate", "intercept_error", "slope_error",
                         "intercept_tval", "slope_tval","intercept_pval", "slope_pval","samplename")
        tmpdf
    })



theta_fitdf = Reduce(function(...) merge(..., all=T), tmplist)
rownames(theta_fitdf) <- theta_fitdf$samplename

    #fitdf <- as.data.frame(t(as.data.frame(fitlist)))
    #colnames(fitdf) <- c("intercept_estimate", "slope_estimate", "intercept_error", "slope_error",
    #                     "intercept_tval", "slope_tval","intercept_pval", "slope_pval")
    for (n in colnames(theta_fitdf)[1:8]){
        theta_fitdf[,n] <- as.numeric(as.character(theta_fitdf[,n]))
    }

theta_fitdf



theta_plot <- ggplot(theta_df, aes(x=yield,y=theta)) + theme_bw(base_size = 18, base_family = "arial") +
      #geom_abline(data=theta_fitdf[theta_fitdf$slope_pval <= 0.05,], 
      #            aes(slope=slope_estimate, intercept=intercept_estimate), alpha=0.5, color="blue") +
      geom_abline(data=theta_fitdf, 
                  aes(slope=slope_estimate, intercept=intercept_estimate, color=samplename), alpha=0.8, size=2) + 

      geom_point(size=5, alpha=1.0) +
      geom_point(size=2, color="LightBlue", alpha=0.8) + 

      scale_x_log10(breaks=xbreaks, labels=xlabs) + 
      scale_y_log10(breaks=ybreaks, labels=ylabs) +  
      #scale_color_manual(values = color_list) +
      xlab(expression("Growth yield : " ~ bgroup("(", over(-Delta*G[cat],Delta*G[ana]+Delta*G[dis]),")"))) +
      ylab(expression(theta ~ "(metacommunity diversity)"))  +
      
      ggtitle("Metacommunity diversity\nacross biomes\n\"Simple\" metabolisms") 
     

theta_plot

theta_fitdf

ggsave("theta_diversity_vs_yield_by_biome_simple_metabolisms.png", plot=theta_plot, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)

ggsave("theta_diversity_vs_yield_biome_colored.png", plot=theta_plot, dpi = 600, 
           units="mm", width = 150*2, height = 100*2, limitsize = FALSE)
#write.csv(theta_fitdf, file="./output_theta/fits_biomlevel.csv", quote=F)

theta_fitdf

#theta_df <- merge(theta_df, biomdf, by.x="dirname", by.y="samplename")
f1 <- lm(log10(theta) ~ log10(yield), data=theta_df)
f3 <- lm(log10(theta) ~ log10(yield)/biom  + biom , data=theta_df)
anova(f1,f3)

theta_df[theta_df$metabolism == "AerobicResp",]

theta_df[theta_df$dirname == "TARA_ERP001736",]

simples  = subset_theta[!(subset_theta$Name %in% c("Aerobic Heterotroph", "Fermenting Heterotroph", "Sulphate reducer")),]


ggplot(data=simples,mapping = aes(x=yield, y=theta) ) + geom_point() + scale_x_log10() + scale_y_log10()

ps <- PSobjs_full[["biom_files//TARA_ERP001736.biom.json"]]

otutab = as.data.frame(otu_table(ps))


relotutab <- apply(otutab, 2, function(x) x/sum(x))

relotutab

mean_abundances = apply(relotutab, 1, mean)

detection_freq = apply(relotutab > 0, 1, function(x) sum(x)/length(x))

piplot_df = data.frame(meanabund = mean_abundances, freq = detection_freq)

piplot_df = piplot_df[piplot_df$meanabund > 1e-6, ]

ggplot(mapping=aes(x=meanabund, y=freq), data=piplot_df) + geom_point(alpha=0.1) + scale_x_log10() + scale_y_log10()

df <- as.data.frame(otu_table(PSobjs_full[["biom_files//TARA_ERP001736.biom.json"]]))

tab = df[rowSums(df)>1,]

obs_freq = apply(tab,1,function(x) sum(x>0) / length(x) * 100) 

mean_abund = apply(apply(tab,2,function(x) x/sum(x)),1,mean) * 100

pidf <- data.frame("of" = obs_freq, "ab"= mean_abund)

(p1<-ggplot(pidf, aes(x=ab, y=of)) + geom_point(alpha=0.3) + scale_x_continuous(limits = c(0,0.025)) + 
theme_bw() + 
    scale_x_log10(limits=c(1e-5,120), breaks=c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1,10,100), labels=c("0.00001", "0.0001","0.001", "0.01","0.1","1","10","100")) +
    labs(x="Mean relative abundance (%)", y="Detection frequency (%)") + guides(color=FALSE) +
    ggtitle('TARA Ocean Sampling\nDetection frequency vs mean abundance'))

ggsave("TARA  plot detection freqs.png", p1, dpi = 600)

for (n in names(PSobjs_full)){
    curname=gsub('biom_files\\/\\/(.*)\\.biom\\.json','\\1', n)
    
    df <- as.data.frame(otu_table(PSobjs_full[[n]]))

    tab = df[rowSums(df)>1,]

    obs_freq = apply(tab,1,function(x) sum(x>0) / length(x) * 100) 

    mean_abund = apply(apply(tab,2,function(x) x/sum(x)),1,mean) * 100

    pidf <- data.frame("of" = obs_freq, "ab"= mean_abund)
    (p1<-ggplot(pidf, aes(x=ab, y=of)) + geom_point(alpha=0.3) + scale_x_continuous(limits = c(0,0.025)) + 
theme_bw() + 
    scale_x_log10(limits=c(1e-5,120), breaks=c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1,10,100), labels=c("0.00001", "0.0001","0.001", "0.01","0.1","1","10","100")) +
    labs(x="Mean relative abundance (%)", y="Detection frequency (%)") + guides(color=FALSE) +
    ggtitle(paste(curname,'\nDetection frequency vs mean abundance', sep = ""))
     )
     
    ggsave(paste("PIPLots/",curname, "  plot detection freqs.png", sep=""), p1, dpi = 600)
}

curname=gsub('biom_files\\/\\/(.*)\\.biom\\.json','\\1', n)
