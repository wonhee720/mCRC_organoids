# Rscript for plotting circos
# Changed the colors to fit with the SV colors from the phylogenetic tree - LWH 22.06.02
# Modified for L1 circos

library(stringr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(circlize)

args <- commandArgs(trailingOnly = TRUE)
L1_path <- args[1]
output_path <- args[2]

counter_CHR = c(
                "chr1"="chr10",
                "chr2"="chr12",
                "chr3"="chr13",
                "chr4"="chr15",
                "chr5"="chr18",
                "chr6"="chr19",
                "chr7"="chrX",
                "chr8"="chrY",
                "chr9"="chrY",
                "chr10"="chr1",
                "chr11"="chr1",
                "chr12"="chr2",
                "chr13"="chr3",
                "chr14"="chr3",
                "chr15"="chr4",
                "chr16"="chr5",
                "chr17"="chr5",
                "chr18"="chr6",
                "chr19"="chr6",
                "chr20"="chr6",
                "chr21"="chr7",
                "chr22"="chr7",
                "chrX"="chr7",
                "chrY"="chr8"
    )
    
    half_POS = c(
                "chr1"=124625310,
                "chr2"=121599686,
                "chr3"=99011215,
                "chr4"=95577138,
                "chr5"=90457630,
                "chr6"=85557534,
                "chr7"=79569332,
                "chr8"=73182011,
                "chr9"=70606716,
                "chr10"=67767374,
                "chr11"=67503258,
                "chr12"=66925948,
                "chr13"=57584939,
                "chr14"=53674770,
                "chr15"=51265696,
                "chr16"=45177376,
                "chr17"=40597605,
                "chr18"=39038624,
                "chr19"=29564492,
                "chr20"=31512760,
                "chr21"=24064948,
                "chr22"=25652283,
                "chrX"=77635280,
                "chrY"=29686783
    )
    
    L1_source = c(
                'chr22_29062287'='22q12.1-2',
                'chr14_59220393'='14q23.1',
                'chr12_117814453'='12q24.22',
                'chr12_71020040'='12q15',
                'chr6_13191024'='6p24.1',
                'chr8_73790808'='8q13.3',
                'chr1_119397988'='1p12',
                'chr11_105478991'='11q22.3',
                'chr14_71197788'='14q24.2',
                'chr1_116980834'='1p13.1',
                'chr2_191478707'='2q32.2',
                'chr6_24814920'='6p22.3-3',
                'chrX_11728384'='Xp22.2-1',
                'chr7_89944417'='7q21.13-1',
                'chr6_53552932'='6p12.1',
                'chr13_30218843'='13q12.3',
                'chr12_3611378'='12p13.32',
                'chr9_5491412'='9p24.1',
                'chr2_11139246'='2p25.1',
                'chr7_90331914'='7q21.13-2'
    )

circos <- function(L1_r){

    chrlist=paste0("chr",c(1:22,"X","Y"))
    L1_pal = c(
       '22q12.1-2'= '#d60000',
       '14q23.1'= '#8c3bff', 
       '12q24.22'=  '#018700', 
       '12q15'='#97ff00', 
       '6p24.1'=   '#00acc6', 
       '8q13.3'='#ff7ed1',
       '1p12'='#6b004f', 
       '11q22.3'='#ffa52f', 
       '14q24.2'='#573b00', 
       '1p13.1'='#005659',
       '2q32.2'='#0000dd',
       '6p22.3-3'='#00fdcf',
       'Xp22.2-1'='#a17569',
       '7q21.13-1'='#bcb6ff',
       '6p12.1'='#95b577',
       '13q12.3'='#bf03b8',
       '12p13.32'= '#645474',
       '9p24.1'='#790000',
       '2p25.1'='#0774d8',
       '7q21.13-2'='#fdf490'
    ) 

    circos.clear()
    circos.par("start.degree" = 90, "gap.degree" = c(rep(1,length(chrlist)-1), 8), canvas.xlim=c(-1.1,1) )
    circos.initializeWithIdeogram(plotType = c("ideogram", "labels"))
            ### L1 ###
    
    baseline1=get.cell.meta.data("cell.bottom.radius")
    baseline2=get.cell.meta.data("cell.bottom.radius")
    height1=0.25
    height2=0.25
    print(baseline1)
    L1_r_TD = L1_r[!is.na(L1_r$src_CHR), ] # Dataframe with identified L1 sources
    L1_r_rest = L1_r[is.na(L1_r$src_CHR), ]  # Dataframe without identified L1 sources 


    # Non-TD (L1s with source unidentified) # Draw them as arrows with only arrow head visible
    if(length(L1_r_rest$CHR)){
        for (i in 1:nrow(L1_r_rest)){
            
            h=NULL
            
            start_CHR = counter_CHR[L1_r_rest$CHR[i]]
            start_POS = half_POS[start_CHR]

            circos.link(sector.index1 = start_CHR, 
                        point1 = start_POS, 
                        sector.index2 = L1_r_rest$CHR[i], 
                        point2 = L1_r_rest$POS[i],
                        #col = '#FFFFFF', # Setting this to white makes the arrow color to white as well. 
                        # So we could only set the color as black and line width 0. 
                        # BUT, in case of vector map, line is visible even with width 0. So we need to render this image as PNG or erase the lines manually in illustrator. -LWH
                        h=h,
                        lwd=0, 
                        rou=baseline1,
                        arr.type = "triangle",
                        arr.length = 0.4,
                        arr.width = 0.4,
                        arr.lwd = 3,
                        arr.col = '#000000', # black
                        directional = 1,
                    )
        }
    }

    # TD (Arrows indicating source to target site transduction)
    if(length(L1_r_TD$CHR)){
        for (i in 1:nrow(L1_r_TD)){
            
            L1_r_TD$src_fin <- paste0(L1_r_TD$src_CHR,"_",L1_r_TD$src_POS)
          
            h=NULL
            Col=L1_pal[L1_source[L1_r_TD$src_fin[i]]]

            circos.link(sector.index1 = L1_r_TD$src_CHR[i], 
                        point1 = L1_r_TD$src_POS[i], 
                        sector.index2 = L1_r_TD$CHR[i], 
                        point2 = L1_r_TD$POS[i],
                        col = Col,
                        h=h,
                        lwd=15, 
                        rou=baseline1,
                        arr.type = "triangle",
                        arr.length = 1,
                        arr.width = 1,
                        arr.lwd = 18,
                        arr.col = Col,
                        directional = 1,
                    )
        }
    }
    
    
    
            ### Legends ###
    
   # lgd4=legend("bottomleft", lty=1, lwd=3, legend=names(L1_pal), col=L1_pal, bty='n', title="Source", title.adj=0.3)
    
}

#png(output_path, res=200, width=1800, height=1600)
pdf(output_path)

L1_r = read_tsv(file=L1_path, col_types=cols(CHR = col_character(), POS = col_double(), src_CHR = col_character(), src_POS = col_double(), L1_type = col_character()))
circos(L1_r)

dev.off()