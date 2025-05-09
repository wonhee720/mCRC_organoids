# Rscript for plotting circos

library(stringr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(circlize)

args <- commandArgs(trailingOnly = TRUE)
snv_path <- args[1]
indel_path <- args[2]
cnv_path <- args[3]
sv_path <- args[4]
output_path <- args[5]
sample_id <- args[6]


circos <- function(snv_r, indel_r, cnv_r, sv_r, text){

    chrlist=paste0("chr",c(1:22,"X","Y"))
    snp_pal = c("C>A"="#15A0EC", # blue
               "C>G"="#0D0C1B", # black
               "C>T"="#F23A29", # red
               "T>A"="#A1A1A1", # grey
               "T>C"="#5AB440", # green
               "T>G"="#ffc001") # yellow
    indel_pal = c("del"="#EC651A", "ins"="#0D33A6")
    line_pal = c("Depth ratio" = "red")
    sv_pal = c("5to5"="#F23A29",
                "3to3"="#F2B705",
                "3to5"="#15A0EC",
                "5to3"="#39A649",
                "translocation"="#DD4DF0")
    sv_pal2 = c("deletion"="#15A0EC",
                "duplication"="#39A649",
                "head-to-head inversion"="#F23A29",
                "tail-to-tail inversion"="#F2B705",
                "translocation"="#DD4DF0")
    maxdiff = ceiling(max(snv_r$logDIFF))

    
    if (max(cnv_r$CNt)>10) {
        maxcnv = 10
    } else {
        maxcnv = max(cnv_r$CNt)
    }
    

    circos.clear()
    circos.par("start.degree" = 90, "gap.degree" = c(rep(1,length(chrlist)-1), 8), canvas.xlim=c(-1.1,1) )
    circos.initializeWithIdeogram()
    
    
            ### SNV ###

    circos.track(factors=chrlist, ylim=c(0,maxdiff),track.height = 0.2) # initialize empty track
    for(CHR in chrlist) {
        # grey lines in background
        for(lvl in seq(0,maxdiff,2)){
            circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="gray", sector.index=CHR)
        }
        circos.points(x=snv_r[snv_r$chr==CHR,]$POS, y=snv_r[snv_r$chr==CHR,]$logDIFF, col=snp_pal[snv_r[snv_r$chr==CHR,]$type],
                      sector.index=CHR, pch=20, cex=0.5, track.index=3)
    }
    circos.yaxis(side="left", at=seq(0, maxdiff, 2), sector.index="chr1", labels.cex=0.6, tick.length = 0.1, track.index = 3)
    
    
            ### Indel ###
    
    circos.track(factors=chrlist, ylim=c(0,1),track.height = 0.1) # initialize empty track
    for(CHR in chrlist) {
        # grey lines in background
        for(lvl in seq(0,1,0.5)){
            circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="gray", sector.index=CHR)
        }
        circos.points(x=indel_r[indel_r$chr==CHR,]$POS, y=indel_r[indel_r$chr==CHR,]$vaf, col=indel_pal[indel_r[indel_r$chr==CHR,]$type],
                      sector.index=CHR, pch=20, cex=0.5, track.index = 4)
    }
    circos.yaxis(side="left", at=seq(0, 1, 0.5), sector.index="chr1", labels.cex=0.6, tick.length = 0.1, track.index = 4)
    
    
            ### CNV ###
    
    circos.track(factors=chrlist, ylim=c(0,maxcnv), track.height = 0.2)
    for(CHR in chrlist) {
        for(lvl in seq(0,maxcnv,2)){
            circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="gray", sector.index=CHR)
        }
        circos.segments(x0=cnv_r[cnv_r$chr==CHR,]$start, x1=cnv_r[cnv_r$chr==CHR,]$end, y0=cnv_r[cnv_r$chr==CHR,]$A, y1=cnv_r[cnv_r$chr==CHR,]$A, 
                          sector.index=CHR, lty=1, lwd=1.5, col="red", track.index = 5) 
                            # Width differene between A and B allele for effective visualization
        circos.segments(x0=cnv_r[cnv_r$chr==CHR,]$start, x1=cnv_r[cnv_r$chr==CHR,]$end, y0=cnv_r[cnv_r$chr==CHR,]$B, y1=cnv_r[cnv_r$chr==CHR,]$B, 
                          sector.index=CHR, lty=1, lwd=1, col="blue", track.index = 5)
        }
    circos.yaxis(side="left", at=seq(0, maxcnv, 2), sector.index="chr1", labels.cex=0.6, tick.length = 0.1, track.index = 5)
    
    
            ### SV ###
    
    baseline1=get.cell.meta.data("cell.bottom.radius")
    baseline2=get.cell.meta.data("cell.bottom.radius")
    height1=0.1
    height2=0.1
    if(length(sv_r$CHR1)){
        for (i in 1:nrow(sv_r)){
            if(sv_r$CHR1[i] == sv_r$CHR2[i]){
                h=ifelse(sv_r$CT[i] %in% c("5to5","3to3"), height1, height2)
                Col=sv_pal[sv_r$CT[i]]
            } else {
                h=NULL
                Col=sv_pal["translocation"]
            }
            circos.link(sector.index1 = sv_r$CHR1[i], point1=sv_r$POS1[i], 
                    sector.index2 = sv_r$CHR2[i], point2=sv_r$POS2[i],
                    col = Col,
                    h=h,
                    lwd=0.7,
                    rou=ifelse(sv_r$CT[i] %in% c("5to5","3to3"), baseline1, baseline2))
        }
    }
    
            ### Legends ###
    
    
    lgd1=legend("bottomright", lty=1, lwd=2, legend=names(line_pal), col=line_pal, bty='n', title="Copy Number")
    lgd2=legend("topleft", pch=21, legend=names(snp_pal), pt.bg=snp_pal, bty='n', title="SNV")
    lgd3=legend("topright", pch=21, legend=names(indel_pal), pt.bg=indel_pal, bty='n', title="Indel")
    lgd4=legend("bottomleft", lty=1, lwd=2, legend=names(sv_pal2), col=sv_pal2, bty='n', title="SV", title.adj=0.3)
    
    mtext(side = 3, text, line = -1, adj = 0.5, cex = 2, font = 2)
        
}

png(output_path, res=200, width=1800, height=1600)

snv_r = read_tsv(file=snv_path, col_types=cols(chr = col_character(), POS = col_double(), logDIFF = col_double(), type = col_character(), context_3 = col_character()))
indel_r = read_tsv(file=indel_path, col_types=cols(chr = col_character(), POS = col_double(), REF = col_character(), ALT = col_character(), type = col_character()))
cnv_r = read_tsv(file=cnv_path, col_types=cols(chr = col_character(), start = col_double(), end = col_double(), Mean = col_double()))   
sv_r = read_tsv(file=sv_path, col_types=cols(CHR1 = col_character(), POS1 = col_double(), CHR2 = col_character(), POS2 = col_double(), CT = col_character(), SVTYPE = col_character()))
circos(snv_r,indel_r,cnv_r,sv_r,sample_id)

dev.off()