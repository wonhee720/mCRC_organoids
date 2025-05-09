circos <- function(snv, cnv, sv, text,
                   snp_pal = c("C>A"="#15A0EC","G>T"="#15A0EC", # blue
                               "C>G"="#0D0C1B","G>C"="#0D0C1B", # black
                               "C>T"="#F23A29","G>A"="#F23A29", # red
                               "T>A"="#A1A1A1","A>T"="#A1A1A1", # grey
                               "T>C"="#5AB440","A>G"="#5AB440", # green
                               "T>G"="#F2BBC5","A>C"="#F2BBC5"), # pink
                   line_pal = c("5'-5'"="#F23A29",
                                "3'-3'"="#F2B705",
                                "3'-5'"="#15A0EC",
                                "5'-3'"="#39A649",
                                "translocation"="#DD4DF0"),
                   chrlist=paste0("chr",c(1:22,"X","Y"))){

    circos.clear()
    circos.par("start.degree" = 90, # start from 12 oclock direction
                 cell.padding = c(0.00, 1.00, 0.00, 1.00), # exact ylim will be applied
                 track.margin = c(0.015, 0.015), # increased gap between track
                 gap.degree = c(rep(1,length(chrlist)-1),8))
    circos.initializeWithIdeogram(plotType = NULL, chromosome.index = chrlist)
    


    # 0. Ideogram ------------------------------------------------------------------------------
    circos.genomicIdeogram(track.height = 0.05)
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        chr = get.cell.meta.data("sector.index") %>% str_replace("chr","")
        xcenter = get.cell.meta.data("xcenter")
        ycenter = get.cell.meta.data("ylim")[2]
        circos.text(xcenter, ycenter, chr,adj = c(0.5,-1.0))
        })
    


    # 1-1. SNV (VAF) ------------------------------------------------------------------------------
    maxvaf = round(max(snv$vaf), 1)
    ymax = ifelse(max(snv$vaf) >= maxvaf, maxvaf + 0.1, maxvaf)
    circos.track(factors=chrlist, ylim=c(0,ymax),track.height = 0.10, bg.border="#00000030") # initialize empty track

    for(CHR in chrlist) {
        # grey lines in background
        for(lvl in 1:(ymax*10)){
            lvl <- lvl/10
            circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="grey", sector.index=CHR)
        }
        circos.points(x=snv[snv$chr==CHR,]$start, y=snv[snv$chr==CHR,]$vaf, col=snp_pal[snv[snv$chr==CHR,]$type],
                      sector.index=CHR, pch=20, cex=0.5)
    }
    circos.yaxis(side="left", at=seq(0, ymax, 0.2), sector.index="chr1", labels.cex=0.6, tick.length = 0.1)


    # 1-2. SNV (logDIFF) ------------------------------------------------------------------------------
    tmp = snv[!is.na(snv$logDIFF), ]
    maxdiff = ceiling(max(tmp$logDIFF))
    mindiff = floor(min(tmp$logDIFF))
    circos.track(factors=chrlist, ylim=c(mindiff,maxdiff),track.height = 0.10, bg.border="#00000030") # initialize empty track

    for(CHR in chrlist) {
        # grey lines in background
        for(lvl in 1:(maxdiff/2)){
            lvl <- lvl*2
            circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="grey", sector.index=CHR)
        }
        circos.points(x=tmp[tmp$chr==CHR,]$start, y=tmp[tmp$chr==CHR,]$logDIFF, col=snp_pal[tmp[tmp$chr==CHR,]$type],
                      sector.index=CHR, pch=20, cex=0.5)
    }
    circos.yaxis(side="left", at=seq(mindiff, maxdiff, 4), sector.index="chr1", labels.cex=0.6, tick.length = 0.1)


    # 2. CNV (total, major, minor CN) & SV (1. middle circle) ---------------------------------------
    cnv_ymax = ceiling(max(cnv$CNt))+1
    circos.track(factors=chrlist, ylim=c(0,cnv_ymax), track.height = 0.2, bg.border="#00000030",bg.col="#90efa810")

    if(length(sv$chr1)){
        sv <- sv %>% filter((chr1 %in% chrlist) & (chr2 %in% chrlist))
        for (i in 1:nrow(sv)){
            if(sv$chr1[i] == sv$chr2[i]){
                Col=line_pal[sv$ori[i]]
            } else {
                Col=line_pal["translocation"]
            }
            circos.segments(sector.index = sv$chr1[i],
                        x0 = sv$pos1[i], y0 = 0,
                        x1 = sv$pos1[i], y1 = cnv_ymax,
                        col = Col,
                        lwd=0.7)
            circos.segments(sector.index = sv$chr2[i],
                        x0 = sv$pos2[i], y0 = 0,
                        x1 = sv$pos2[i], y1 = cnv_ymax,
                        col = Col,
                        lwd=0.7)
        }
    }


    for(CHR in chrlist) {
        for(lvl in 1:(cnv_ymax-1)){
            circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="grey70", sector.index=CHR)
        }
        circos.points(x=cnv[cnv$chr==CHR,]$start, y=cnv[cnv$chr==CHR,]$CNt, 
                      sector.index=CHR, pch=20, cex=0.5, col="black", )
    }
    circos.yaxis(side="left", at=seq(0, cnv_ymax, 1), sector.index="chr1", labels.cex=0.6, tick.length = 0.1)



    # 3. SV (2. inner-most circle) -------------------------------------------------------------------------
    baseline1=get.cell.meta.data("cell.bottom.radius")
    baseline2=get.cell.meta.data("cell.bottom.radius")
    height1=0.1
    height2=0.1
    print(sv$chr1)
    if(length(sv$chr1)){
        for (i in 1:nrow(sv)){
            if(sv$chr1[i] == sv$chr2[i]){
                h=ifelse(sv$ori[i] %in% c("5'-5'","3'-3'"), height1, height2)
                Col=line_pal[sv$ori[i]]
            } else {
                h=NULL
                Col=line_pal["translocation"]
            }
            circos.link(sector.index1 = sv$chr1[i], sv$pos1[i], 
                    sector.index2 = sv$chr2[i], sv$pos2[i],
                    col = Col,
                    h=h,
                    lwd=0.7,
                    rou=ifelse(sv$type[i] %in% c("5'-5'","3'-3'"), baseline1, baseline2))
        }
    }



    mtext(side = 3, text, outer = F, line = -1, adj = 0.00, cex = 1.1, font = 4)
    snp_pal = c("C>A"="#15A0EC", # blue
                "C>G"="#0D0C1B", # black
                "C>T"="#F23A29", # red
                "T>A"="#A1A1A1", # grey
                "T>C"="#5AB440", # green
                "T>G"="#F2BBC5") # pink
    line_pal = c("5'-5'"="#F23A29",
                 "3'-3'"="#F2B705",
                 "3'-5'"="#15A0EC",
                 "5'-3'"="#39A649",
                 "translocation"="#DD4DF0")

    lgd1=legend("bottomleft", lty=1, legend=names(line_pal), col=line_pal, bty='n')
    lgd2=legend("bottomright", legend=c("major CN","minor CN"), bty='n',col=c("#EC651A","#0D33A6"), lty=1, lwd=3)
    lgd3=legend(lgd2$rect$left + lgd2$rect$w, lgd2$rect$top - 0.05, xjust = 1, yjust=0, pch=21, legend=names(snp_pal), pt.bg=snp_pal, bty='n')

}
