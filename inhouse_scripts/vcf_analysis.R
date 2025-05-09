library(tidyverse)
library(cowplot)
library(mclust)
setwd("/home/users/kimin/projects/03_Mucosal_Melanoma/")
list.files()

vcf <- read_tsv("02_vcf/MM001-Tumor-Blood_mns.tsv", col_types = cols('CHROM' = 'c'))
vaf <- vcf$TUMOR_AF

summary(vcf['TUMOR_AF'])
mean(vaf)
median(vaf)

vafmode <- function(vec){
  mt <- table(vec)
  names(mt)[mt == max(mt)]
}
vafmode(vaf)

hist(vaf)
plot(density(vaf, kernel = "gaussian", adjust = 1))

pdf(file = "practice.pdf", width = 11)
par(mar = c(5,5,2,2))
plot(density(vaf, kernel = "gaussian", adjust = 1))
lines(c(0.067,0.067), c(0,25), type='l')
dev.off()



conv_v=rep(c("CA","CG","CT","TA","TC","TG"),2)
name_v=c("CA","CG","CT","TA","TC","TG","GT","GC","GA","AT","AG","AC")
names(conv_v)=name_v
context_color=c("#1EBFF0", "#050708", "#E62725", "#CBCACB", "#A1CF64", "#EDC8C5") 



#filtering with defined cutoff
#var_readN-11 > log2(OR+0.001)
file_list <- list.files(pattern = "anv.fi2$")
length(file_list)
#dtlist_woORfi <- list()
for (file in file_list){
  sampleid <- unlist(strsplit(file,'\\.'))[1]
  dt <- read_tsv(file, col_types = cols(`#CHROM`="c"))
  colnames(dt)
  dt <- dt %>% separate("pcawgN36_DP;DP_N;Ref;RefN;Var;VarN;Error_vaf", c("pcawgN36_DP","pcawgN36_DP_N","pcawgN36_Ref","pcawgN36_RefN","pcawgN36_Var","pcawgN36_VarN","pcawgN36_Error_vaf"), sep=';', remove=F) %>%
    separate(`snuN30_DP;DP_N;Ref;RefN;Var;VarN;Error_vaf`, c('snuN30_DP','snuN30_DP_N','snuN30_Ref','snuN30_RefN','snuN30_Var','snuN30_VarN','snuN30_Error_vaf'), sep=';', remove=F) %>%
    separate(`bgiN24_DP;DP_N;Ref;RefN;Var;VarN;Error_vaf`, c('bgiN24_DP','bgiN24_DP_N','bgiN24_Ref','bgiN24_RefN','bgiN24_Var','bgiN24_VarN','bgiN24_Error_vaf'), sep= ';', remove=F)
  dt <- dt %>% mutate_at(c("pcawgN36_DP","pcawgN36_DP_N","pcawgN36_Ref","pcawgN36_RefN","pcawgN36_Var","pcawgN36_VarN","pcawgN36_Error_vaf"), funs(as.numeric(.))) %>%
    mutate_at(c('snuN30_DP','snuN30_DP_N','snuN30_Ref','snuN30_RefN','snuN30_Var','snuN30_VarN','snuN30_Error_vaf'), funs(as.numeric(.))) %>%
    mutate_at(c('bgiN24_DP','bgiN24_DP_N','bgiN24_Ref','bgiN24_RefN','bgiN24_Var','bgiN24_VarN','bgiN24_Error_vaf'), funs(as.numeric(.)))
  #  dt <- dt %>% mutate(snuN30_var_R = snuN30_Var/snuN30_VarN, pcawgN36_var_R = pcawgN36_Var/pcawgN36_VarN, bgiN24_var_R = bgiN24_Var/bgiN24_VarN) %>%
  #    mutate(snuN30_Var_new = ifelse(snuN30_Var > 0 & snuN30_var_R >= 4, 0, snuN30_Var), pcawgN36_Var_new = ifelse(pcawgN36_Var > 0 & pcawgN36_var_R >=4, 0, pcawgN36_Var), bgiN24_Var_new = ifelse(bgiN24_Var > 0 & bgiN24_var_R >=4, 0, bgiN24_Var))
  dt <- dt %>% mutate(total_ref =snuN30_Ref + pcawgN36_Ref + bgiN24_Ref, total_var=snuN30_Var + pcawgN36_Var + bgiN24_Var)
  dt <- dt %>% select(-c("pcawgN36_DP","pcawgN36_DP_N","pcawgN36_Ref","pcawgN36_RefN","pcawgN36_Var","pcawgN36_VarN","pcawgN36_Error_vaf")) %>% select(-c('snuN30_DP','snuN30_DP_N','snuN30_Ref','snuN30_RefN','snuN30_Var','snuN30_VarN','snuN30_Error_vaf')) %>%
    select(-c('bgiN24_DP','bgiN24_DP_N','bgiN24_Ref','bgiN24_RefN','bgiN24_Var','bgiN24_VarN','bgiN24_Error_vaf'))
  dt[is.na(dt)]<-0
  dt$pvalue <- dt[c("ref_readN","var_readN","total_ref","total_var")] %>% 
    apply(1, function(x){tbl <- matrix(as.numeric(x[1:4]), nrow=2); fisher.test(tbl)$p.value})
  dt$OR <- dt[c("ref_readN","var_readN","total_ref","total_var")] %>% 
    apply(1, function(x){tbl <- matrix(as.numeric(x[1:4]), nrow=2); fisher.test(tbl)$estimate})
  #dt <- dt %>% mutate(context = as.character(conv_v[paste0(REF,ALT)]))
  dt <- dt %>% mutate(context = as.character(conv_v[paste0(REF,ALT)]))
  #dtlist_woORfi[[sub("\\..*","",file)]] <- dt
  #dt <- dt[is.na(dt$POS) == F,]
  #png('tmp.png',width=16, height=9, res=300, units="in")
  #g <- ggplot(dt_1000_1, aes(x=var_readN, y=log2(OR+0.001)))+
  #  geom_jitter(aes(color=context),alpha=0.5, width=0.1, height=0)+
  #  geom_segment(aes(x=0, xend=10, y=-11, yend=-2), linetype = "dotted", color = "red", alpha=0.5)+
  #  geom_segment(aes(x=0, xend=30, y=-5, yend=-5), linetype = "dotted", color = "red", alpha=0.5)+
  #  scale_x_continuous(limits=c(0,30))+
  #  scale_y_continuous(breaks=seq(-10,-2,1), limits=c(-11.5,-1.5))+
  #  ggtitle("dtf1")
  #print(g)
  #dt %>% filter(var_readN -11 > log2(OR + 0.001)& log2(OR + 0.001) < -5) %>% write_tsv(paste0(file,'.np_fi'))
  #dt %>% filter( log2(OR + 0.001) < -6) %>% write_tsv(paste0(file,'.np_fi2'))
}
########


#########

for (file in file_list){
  sampleid <- unlist(strsplit(file,'\\.'))[1]
  dt <- read_tsv(file, col_types = cols(`#CHROM`="c"))
  colnames(dt)
  dt <- dt %>% separate("pcawgN36_DP;DP_N;Ref;RefN;Var;VarN;Error_vaf", c("pcawgN36_DP","pcawgN36_DP_N","pcawgN36_Ref","pcawgN36_RefN","pcawgN36_Var","pcawgN36_VarN","pcawgN36_Error_vaf"), sep=';', remove=F) %>%
    separate(`snuN30_DP;DP_N;Ref;RefN;Var;VarN;Error_vaf`, c('snuN30_DP','snuN30_DP_N','snuN30_Ref','snuN30_RefN','snuN30_Var','snuN30_VarN','snuN30_Error_vaf'), sep=';', remove=F) %>%
    separate(`bgiN24_DP;DP_N;Ref;RefN;Var;VarN;Error_vaf`, c('bgiN24_DP','bgiN24_DP_N','bgiN24_Ref','bgiN24_RefN','bgiN24_Var','bgiN24_VarN','bgiN24_Error_vaf'), sep= ';', remove=F) %>%
    separate(`pcnslN11_DP;DP_N;Ref;RefN;Var;VarN;Error_vaf`, c("pcnslN11_DP","pcnslN11_DP_N","pcnslN11_Ref","pcnslN11_RefN","pcnslN11_Var","pcnslN11_VarN","pcnslN11_Error_vaf"), sep=";", remove=F)
  dt <- dt %>% mutate_at(c("pcawgN36_DP","pcawgN36_DP_N","pcawgN36_Ref","pcawgN36_RefN","pcawgN36_Var","pcawgN36_VarN","pcawgN36_Error_vaf"), funs(as.numeric(.))) %>%
    mutate_at(c('snuN30_DP','snuN30_DP_N','snuN30_Ref','snuN30_RefN','snuN30_Var','snuN30_VarN','snuN30_Error_vaf'), funs(as.numeric(.))) %>%
    mutate_at(c('bgiN24_DP','bgiN24_DP_N','bgiN24_Ref','bgiN24_RefN','bgiN24_Var','bgiN24_VarN','bgiN24_Error_vaf'), funs(as.numeric(.))) %>%
    mutate_at(c("pcnslN11_DP","pcnslN11_DP_N","pcnslN11_Ref","pcnslN11_RefN","pcnslN11_Var","pcnslN11_VarN","pcnslN11_Error_vaf"), funs(as.numeric(.))) 
  #  dt <- dt %>% mutate(snuN30_var_R = snuN30_Var/snuN30_VarN, pcawgN36_var_R = pcawgN36_Var/pcawgN36_VarN, bgiN24_var_R = bgiN24_Var/bgiN24_VarN) %>%
  #    mutate(snuN30_Var_new = ifelse(snuN30_Var > 0 & snuN30_var_R >= 4, 0, snuN30_Var), pcawgN36_Var_new = ifelse(pcawgN36_Var > 0 & pcawgN36_var_R >=4, 0, pcawgN36_Var), bgiN24_Var_new = ifelse(bgiN24_Var > 0 & bgiN24_var_R >=4, 0, bgiN24_Var))
  dt <- dt %>% mutate(total_ref =snuN30_Ref + pcawgN36_Ref + bgiN24_Ref + pcnslN11_Ref, total_var=snuN30_Var + pcawgN36_Var + bgiN24_Var + pcnslN11_Var)
  dt <- dt %>% select(-c("pcawgN36_DP","pcawgN36_DP_N","pcawgN36_Ref","pcawgN36_RefN","pcawgN36_Var","pcawgN36_VarN","pcawgN36_Error_vaf")) %>% select(-c('snuN30_DP','snuN30_DP_N','snuN30_Ref','snuN30_RefN','snuN30_Var','snuN30_VarN','snuN30_Error_vaf')) %>%
    select(-c('bgiN24_DP','bgiN24_DP_N','bgiN24_Ref','bgiN24_RefN','bgiN24_Var','bgiN24_VarN','bgiN24_Error_vaf')) %>% select(-c("pcnslN11_DP","pcnslN11_DP_N","pcnslN11_Ref","pcnslN11_RefN","pcnslN11_Var","pcnslN11_VarN","pcnslN11_Error_vaf"))
  dt[is.na(dt)]<-0
  dt$pvalue <- dt[c("ref_readN","var_readN","total_ref","total_var")] %>% 
    apply(1, function(x){tbl <- matrix(as.numeric(x[1:4]), nrow=2); fisher.test(tbl)$p.value})
  dt$OR <- dt[c("ref_readN","var_readN","total_ref","total_var")] %>% 
    apply(1, function(x){tbl <- matrix(as.numeric(x[1:4]), nrow=2); fisher.test(tbl)$estimate})
  #dt <- dt %>% mutate(context = as.character(conv_v[paste0(REF,ALT)]))
  dt <- dt %>% mutate(context = as.character(conv_v[paste0(REF,ALT)]))
  dt <- dt[is.na(dt$POS) == F,]
  assign(paste0("dt_", sampleid), dt)
  #dt %>% filter(var_readN -11 > log2(OR + 0.001)& log2(OR + 0.001) < -5) %>% write_tsv(paste0(file,'.np_fi'))
  #dt %>% filter( log2(OR + 0.001) < -6) %>% write_tsv(paste0(file,'.np_fi2'))
}

##############################################################

# Loading vcf files with python filter only with additional annotations
dtlist_woORfi <- list()
for (i in idlist) {
  dt1 <- read_tsv(paste0(i, ".snv.mutect2_strelka2_union.vcf.readinfo.readc.rasm.snuN30.pcawgN36.bgiN24.pcnslN11.anv.fi2.ORhi"), 
                  col_types=cols(`#CHROM`="c", "context"="f", "cosmic_count"="f", "gmcluster"="f", "group_by_line"="f", "gmcluster2"="f", "logOR4cutoff"="f"))
  dt2 <- read_tsv(paste0(i, ".snv.mutect2_strelka2_union.vcf.readinfo.readc.rasm.snuN30.pcawgN36.bgiN24.pcnslN11.anv.fi2.ORlo"), 
                  col_types=cols(`#CHROM`="c", "context"="f", "cosmic_count"="f", "gmcluster"="f", "group_by_line"="f", "gmcluster2"="f", "logOR4cutoff"="f"))
  dtlist_woORfi[[i]] <- rbind(dt1, dt2)
  rm(dt1, dt2)
}


#############################################################

# Adding logOR column
for (i in idlist) {
  dtlist_woORfi[[i]] <- dtlist_woORfi[[i]] %>% mutate(logOR=log2(OR+0.001))
}
#####

# kmeans clustering - failed
for (i in sampleid_list) {
  tmp <- mutate( get(paste0("dt_",i)), logOR=log2(OR+0.001) )
  tmp$cluster <- kmeans(select(tmp, c(var_readN, logOR)), 2)$cluster
  assign(paste0("dt_", i), tmp)
  rm(tmp)
}
#####

# segmentation with arbitrary lines ( y = +/- 2/30 * x - 6 )
dt_LY05_C <- mutate (dt_LY05_C, group = if_else( logOR > (-2/30)*var_readN - 6 & logOR < (2/30)*var_readN - 6,
                                                 "mid",
                                                 if_else( logOR >= (2/30)*var_readN - 6 , "top", "bot") ))
#####

# smoothScatter plots
for (i in names(dtlist_woORfi)){
  par(mfrow=c(2,3), mar=c(2,2,2,2))
  for (dn in c('CA', 'CG', 'CT', 'TA', 'TC', 'TG')) {
    dtlist_woORfi[[i]] %>%
      na.omit() %>%
      filter(context == dn) %>%
      select(var_readN, logOR) %>% smoothScatter(main=paste(i, dn), xlim=c(0,40), ylim=c(-11,0), xlab="", ylab="")
  }
}


# PAM using Gower metric. It was failed with large samples (>50000) due to long vector problem(LY01_C) or insufficient memory? problem(LY04_C).
library(cluster)
for (i in dtid_list[-c(1:3)]) {
  p <- get(i) %>%
    select(var_readN, logOR, context) %>%
    daisy(metric="gower", weights=c(2,2,1)) %>%
    pam(diss=T, 3)
  assign(i, mutate(get(i), pamcluster_2_2_1=as.factor(p$clustering)))
  #attributes(get(i))[["silhouette_pam_2_2_1"]] <- p$silinfo$avg.width  # This was failed due to error 'get<- function is not found'. Using get function within attributes function raises error.
  rm(p)
}
###################


# Generation of cosmic_count column. cosmic86_coding column was parsed and all occurrences were added.
myf <- function(x){sum( as.numeric(str_replace(unlist(strsplit(x[2],",")),"\\(.*\\)","")) )}
for (i in idlist){
  dtlist_woORfi[[i]] <- 
    dtlist_woORfi[[i]] %>%
    mutate(cosmic_count=if_else(cosmic86_coding==".", as.numeric(0), unlist(lapply(strsplit(cosmic86_coding, "OCCURENCE="), myf))))
}
####################


# Gaussian mixture model based clustering -- It was not able to cluster LY01_C with subsetting, so LY01_C was operated on without subset.
# All other samples were clustered with initial subsetting during hierarchical clustering.
for (i in dtid_list[-1]){
  i <- "dt_LY01_C"
  tmp <- get(i) %>%
    mutate(rownames=rownames(get(i)), measure=logOR-var_readN) %>%
    select(rownames, measure) %>%
    arrange(measure)
  n <- as.numeric(tmp[c(1:100,(dim(tmp)[1]-99):dim(tmp)[1]),]$rownames)
  set.seed(123)
  mod <- Mclust(select(get(i), var_readN, logOR), G=2, modelNames="VVV", initialization=list(subset=n))
  assign(i, mutate(get(i), gmcluster2=as.factor(mod$classification)))
  rm(mod, tmp, n)
}
##########################



### ggplot making code
for (i in idlist){
  i="LY02_S1"
  cosset <- as.character(unique(dtlist_woORfi[[i]]$cosmic_count))
  cosset <- cosset[order(nchar(cosset), cosset)]
  cosset <- cosset[-1]
  colfunc <- colorRampPalette(c("white", "blue"))
  coscol <- colfunc(length(cosset)+1)[-1]
  names(coscol) <- cosset
  colorset <- c('1'=hcl(0, 100, 60), '2'=hcl(120, 100, 60), '3'=hcl(240, 100, 60), 
                'CA'=context_color[1],'CG'=context_color[2],'CT'=context_color[3],'TA'=context_color[4],'TC'=context_color[5],'TG'="red",
                coscol)
  
  print(
    dtlist_woORfi[[i]] %>% #filter(`#CHROM` %in% c(1,3,9,19)) %>%
      ggplot(aes(x=var_readN, y=logOR)) +
      geom_jitter(aes(color=context), alpha=0.3, width=0.1, height=0) +
      geom_point(aes(color=cosmic_count), data=filter(dtlist_woORfi[[i]], cosmic_count!="coscnt_0"), alpha=1) +
      geom_text(aes(label=str_replace(cosmic_count, "coscnt_", "")), check_overlap=F, data=filter(dtlist_woORfi[[i]], cosmic_count!="coscnt_0"), hjust=0, vjust=0) +
      #geom_segment(aes(x=0, xend=15, y=-7, yend=-3), linetype = "dotted", color = "red", alpha=0.5)+
      #geom_segment(aes(x=0, xend=30, y=-7, yend=-7), linetype = "dotted", color = "red", alpha=0.5)+
      geom_hline(aes(yintercept=-4)) +
      #geom_abline(slope=-2/30, intercept=-6) +
      #geom_abline(slope=2/30, intercept=-6) +
      scale_x_continuous(limits=c(0,50)) +
      #scale_y_continuous(breaks=seq(-10,-2,1), limits=c(-10.5,-1.5)) +
      labs(title=paste(i, "post-filter", "all chrs"))+
      #scale_color_manual(values=c("clonal"=hcl(alpha=0), "early"=rainbow(3)[1], "late/minor"=rainbow(3)[2], "subclonal"=rainbow(3)[3]), breaks=c("clonal", "early", "late/minor", "subclonal"))
      scale_color_manual(values=colorset)+
      theme_bw()+theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line=element_line())
    #, breaks=c('1', '2', names(coscol)), labels=c('1', '2', names(coscol)))
    #theme(legend.position="none")
  )
}
#####
dtlist_wo

which(is.na(dtlist[["LY10_C"]]$timing))
dtlist[["LY10_C"]][c(2271,2272,18203,18204),80:84]

### grid ggplot 
for (i in idlist){
  cosset <- as.character(unique(dtlist_woORfi[[i]]$cosmic_count))
  cosset <- cosset[order(nchar(cosset), cosset)]
  cosset <- cosset[-1]
  colfunc <- colorRampPalette(c("white", "blue"))
  coscol <- colfunc(length(cosset)+1)[-1]
  names(coscol) <- cosset
  colorset <- c('1'=hcl(0, 100, 60), '2'=hcl(120, 100, 60), '3'=hcl(240, 100, 60), 
                'CA'=context_color[1],'CG'=context_color[2],'CT'=context_color[3],'TA'=context_color[4],'TC'=context_color[5],'TG'="red",
                coscol)
  assign(
    paste0("g_",i),
    dtlist_woORfi[[i]] %>%
      ggplot(aes(x=var_readN, y=logOR)) +
      geom_jitter(aes(color=context), alpha=0.1, width=0.1, height=0) +
      geom_point(aes(color=cosmic_count), data=filter(dtlist_woORfi[[i]], cosmic_count!="coscnt_0"), alpha=1) +
      geom_text(aes(label=str_replace(cosmic_count, "coscnt_", "")), check_overlap=F, data=filter(dtlist_woORfi[[i]], cosmic_count!="coscnt_0"), hjust=0, vjust=0) +
      #geom_segment(aes(x=0, xend=15, y=-7, yend=-3), linetype = "dotted", color = "red", alpha=0.5)+
      #geom_segment(aes(x=0, xend=30, y=-7, yend=-7), linetype = "dotted", color = "red", alpha=0.5)+
      geom_hline(aes(yintercept=-4)) +
      #scale_x_continuous(limits=c(0,40)) +
      #scale_y_continuous(breaks=seq(-10,-2,1), limits=c(-10.5,-1.5)) +
      labs(title=i) +
      scale_color_manual(values=colorset, name="", breaks=c('1', '2', names(coscol)), 
                         labels=c('1', '2', names(coscol), guides=NULL)))+
    theme(legend.position="none")
}

i="LY04_C"

plot_grid(plotlist=mget(ls(pattern="g_LY*")[1:6]), nrow=2, ncol=3, align="h")
plot_grid(plotlist=mget(ls(pattern="g_LY*")[7:12]), nrow=2, ncol=3, align="h")
plot_grid(plotlist=mget(ls(pattern="g_LY*")[13:14]), nrow=1, ncol=2, align="h")
rm(list=ls(pattern="g_LY*"))
###


### getting slope and intercept of classification line ==> represented in group_by_line column
x <- as.numeric(NULL)
y <- as.numeric(NULL)
for (i in slopelist){
  tmp <- get(i) %>%
    filter(gmcluster2==2) %>%
    mutate(measure=var_readN+logOR) %>%
    filter(measure==max(measure) | measure==min(measure))
  x <- c(x, tmp[["var_readN"]])
  y <- c(y, tmp[["logOR"]])
  rm(tmp)
}
lm(y~x)$coefficients
#(Intercept)           x 
#-10.734470     0.505705 
###





# Appending corrected var_readN
for (i in idlist){
  dtlist[[i]] <- dtlist[[i]] %>% mutate(var_readN_corrected=var_readN/(filter(dt_bysample, sampleid==i)$vafpur))
  dtlist_woORfi[[i]] <- dtlist_woORfi[[i]] %>% mutate(var_readN_corrected=var_readN/(filter(dt_bysample, sampleid==i)$vafpur))
}

### Timing signature plotting
for (i in idlist[-6]){
  cosset <- as.character(unique(dtlist[[i]]$cosmic_count))
  cosset <- cosset[order(nchar(cosset), cosset)]
  cosset <- cosset[-1]
  colfunc <- colorRampPalette(c("white", "blue"))
  coscol <- colfunc(length(cosset)+1)[-1]
  names(coscol) <- cosset
  colorset <- c('1'=hcl(0, 100, 60), '2'=hcl(120, 100, 60), '3'=hcl(240, 100, 60), 
                'CA'=context_color[1],'CG'=context_color[2],'CT'=context_color[3],'TA'=context_color[4],'TC'=context_color[5],'TG'=context_color[6],
                coscol)
  print(
    dtlist[[i]] %>% filter(timing!="clonal") %>%
      ggplot(aes(x=var_readN, y=logOR)) +
      geom_jitter(aes(color=timing), alpha=0.1, width=0.1, height=0) +
      #geom_point(aes(color=cosmic_count), data=filter(dtlist_woORfi[[i]], cosmic_count!="coscnt_0"), alpha=1) +
      #geom_text(aes(label=str_replace(cosmic_count, "coscnt_", "")), check_overlap=F, data=filter(dtlist_woORfi[[i]], cosmic_count!="coscnt_0"), hjust=0, vjust=0) +
      #geom_segment(aes(x=0, xend=15, y=-7, yend=-3), linetype = "dotted", color = "red", alpha=0.5)+
      #geom_segment(aes(x=0, xend=30, y=-7, yend=-7), linetype = "dotted", color = "red", alpha=0.5)+
      geom_hline(aes(yintercept=-4)) +
      scale_x_continuous(limits=c(0,50)) +
      #scale_y_continuous(breaks=seq(-10,-2,1), limits=c(-10.5,-1.5)) +
      labs(title=paste(i, "timing"))
    #scale_color_manual(values=colorset, name="", breaks=c('1', '2', names(coscol)), labels=c('1', '2', names(coscol)))
    #theme(legend.position="none")
  )
}
#####


### var_readN logOR signature plotting - with var_readN_corrected
for (i in idlist){
  cosset <- as.character(unique(dtlist_woORfi[[i]]$cosmic_count))
  cosset <- cosset[order(nchar(cosset), cosset)]
  cosset <- cosset[-1]
  colfunc <- colorRampPalette(c("white", "blue"))
  coscol <- colfunc(length(cosset)+1)[-1]
  names(coscol) <- cosset
  colorset <- c('1'=hcl(0, 100, 60), '2'=hcl(120, 100, 60), '3'=hcl(240, 100, 60), 
                'CA'=context_color[1],'CG'=context_color[2],'CT'=context_color[3],'TA'=context_color[4],'TC'=context_color[5],'TG'=context_color[6],
                coscol)
  print(
    dtlist_woORfi[[i]] %>% 
      ggplot(aes(x=var_readN_corrected, y=logOR)) +
      geom_jitter(aes(color=context), alpha=0.1, width=0.1, height=0) +
      geom_point(aes(color=cosmic_count), data=filter(dtlist_woORfi[[i]], cosmic_count!="coscnt_0"), alpha=1) +
      geom_text(aes(label=str_replace(cosmic_count, "coscnt_", "")), check_overlap=F, data=filter(dtlist_woORfi[[i]], cosmic_count!="coscnt_0"), hjust=0, vjust=0) +
      #geom_segment(aes(x=0, xend=15, y=-7, yend=-3), linetype = "dotted", color = "red", alpha=0.5)+
      #geom_segment(aes(x=0, xend=30, y=-7, yend=-7), linetype = "dotted", color = "red", alpha=0.5)+
      #geom_hline(aes(yintercept=-4)) +
      scale_x_continuous(limits=c(0,150)) +
      scale_y_continuous(breaks=seq(-10,-2,1), limits=c(-10.5,0)) +
      labs(title=paste(i, "var_readN_corrected"), color="context") +
      scale_color_manual(values=colorset, breaks=c('1', '2', names(coscol)), labels=c('1', '2', names(coscol)))
    #theme(legend.position="none")
  )
}
#####




#OR 높은것 일부, 낮은것 일부 뽑아서 IGV feature, context, signature 등을 비교해 보기
#OR cutoff를 잡는 것이 목적
#OR 구할때 pcnsl 정보 더해서 total_var 만들어서 plot 다시 그리기
