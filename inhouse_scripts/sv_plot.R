### Structural Variation plot with absCN, VAFs
# 2019.05.29 cjyoon
# 2019.07.01 cjyoon added mutated copy number calculation
# 2019.07.28 cjyoon added intermutation distance
# 2019.08.18 cjyoon added kataegis mark
# 2019.10.08 cjyoon argument parsing for snakefile use
# 2019.10.25 cjyoon ideogram added
suppressMessages(library(grid))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(grideogram))
option_list = list(
  make_option(c("--snv_vaf"), type="character", default=NULL,
              help="snv_vaf text file", metavar="character"),
  make_option(c("-s", "--sample_name"), type="character", default=".",
              help="sample name", metavar="character"), 
  make_option(c("--sv"), type="character", default=".",
              help="SV text file", metavar="character"),
  make_option(c("--abscn"), type="character", default=".",
              help="abscn file", metavar="character"), 
  make_option(c("-o", "--output_dir"), type="character", default=".",
              help="output folder path [default= %default]", metavar="character"),
  make_option(c("-r", "--reference"), type="character", default=".",
              help="reference fasta file, needs to be indexed", metavar="character")
  
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

reference = opt$reference 
sample_name = opt$sample_name
snv_vaf_path = opt$snv_vaf
abscn_path = opt$abscn
sv_path = opt$sv
output_dir = opt$output_dir

parse_chrom_size <- function(reference_path){
  fasta_index = paste0(reference_path, '.fai')
  reference=read_delim(fasta_index, delim='\t', col_names=c('chromosome', 'size', 'a', 'b', 'c'))
  structure(reference$size, names=reference$chromosome)
}
# chromosome sizes
chrom_sizes = parse_chrom_size(opt$reference)
chrom_sizes = chrom_sizes[str_detect(names(chrom_sizes), pattern = "^[0-9XY]*$")] # only get autosomal + X Y chromosomes, ignore small contigs
# Simplify to x.xx MB
simplify_coordinate <- function(x) {
  digits = floor(log10(x)) + 1
  # print(digits)
  if (digits >= 6) {
    return(paste0(round(x / 10 ** 6, 2), 'mb'))
  }
  else if (digits >= 3) {
    return(paste0(round(x / 10 ** 3, 2), 'kb'))
  }
  else
    (return(x))
}

simplify_coordinate_vec = Vectorize(simplify_coordinate)

# Function to draw SVs
sv_curve <-
  function(bp1_chrom,
           bp1,
           bp2_chrom,
           bp2,
           start,
           range_len,
           svtype,
           direction) {

    sv_color = c(
      'BND' = 'black', 
      'DEL' = 'blue', 
      'DUP' = 'red', 
      'INV' = '#E3D200', # dyellow in circos
      'INS' = 'orange', 
      'TRA' ='black')

    if (svtype !="BND") {
      plot_start = (bp1 - start) / range_len
      plot_end = (bp2 - start) / range_len
      # cat(plot_start,  plot_end)
      x = 0:100 / 100
      if (svtype == 'DEL') {
        y = -4 * x * (x - 1)
        pushViewport(viewport(
          x = plot_start,
          y = 2 / 3,
          width = (plot_end - plot_start),
          height = 0.15,
          just = c('left', 'bottom'),
          clip = 'off'
        ))
        grid.lines(
          x = c(0, 0),
          y = c(0,-10),
          gp = gpar(col = sv_color[['DEL']], lwd = 0.1)
        )
        grid.lines(
          x = c(1, 1),
          y = c(0,-10),
          gp = gpar(col = sv_color[['DEL']], lwd = 0.1)
        )
        grid.lines(x = x,
                   y = y,
                   gp = gpar(col = sv_color[['DEL']], alpha = 0.5))
        upViewport()
      } else if (svtype == 'DUP') {
        y = 4 * (x - 0.5) ** 2
        pushViewport(viewport(
          x = plot_start,
          y = 2 / 3,
          width = (plot_end - plot_start),
          height = 0.15,
          just = c('left', 'top'),
          clip = 'off'
        ))
        grid.lines(x = x,
                   y = y,
                   gp = gpar(col = sv_color[['DUP']],  alpha = 0.5))
        grid.lines(
          x = c(0, 0),
          y = c(1,-9),
          gp = gpar(col = sv_color[['DUP']], lwd = 0.1)
        )
        grid.lines(
          x = c(1, 1),
          y = c(1,-9),
          gp = gpar(col = sv_color[['DUP']], lwd = 0.1)
        )

        upViewport()
        
      } else if (svtype == 'INV' && direction == '5to5') {
        y = -4 * x * (x - 1)
        pushViewport(viewport(
          x = plot_start,
          y = 1 / 3,
          width = (plot_end - plot_start),
          height = 0.15,
          just = c('left', 'bottom'),
          clip = 'off'
        ))
        grid.lines(
          x = c(0, 0),
          y = c(0,-8),
          gp = gpar(col = sv_color[['INV']], lwd = 0.1)
        )
        grid.lines(
          x = c(1, 1),
          y = c(0,-8),
          gp = gpar(col =  sv_color[['INV']], lwd = 0.1)
        )
        grid.lines(x = x,
                   y = y,
                   gp = gpar(col =  sv_color[['INV']], alpha = 0.5))
        upViewport()
        
      } else if (svtype == 'INV' && direction == '3to3') {
        y = 4 * (x - 0.5) ** 2
        # grid.rect()
        pushViewport(viewport(
          x = plot_start,
          y = 1 / 3,
          width = (plot_end - plot_start),
          height = 0.15,
          just = c('left', 'top'),
          clip = 'off'
        ))
        grid.lines(
          x = c(0, 0),
          y = c(1,-7),
          gp = gpar(col =  sv_color[['INV']], lwd = 0.1)
        )
        grid.lines(
          x = c(1, 1),
          y = c(1,-7),
          gp = gpar(col =  sv_color[['INV']], lwd = 0.1)
        )
        grid.lines(
          x = x,
          y = y,
          gp = gpar(col =  sv_color[['INV']], alpha = 0.5)
        )
        upViewport()
      }
    }
    else{
      # If SV is a translocation
      if (bp1_chrom == chrom) {
        plot_start = (bp1 - start) / range_len
        plot_end = plot_start + 0.1
        mate_chrom = bp2_chrom
        mate_pos = bp2
        
      } else if (bp2_chrom == chrom) {
        plot_start = (bp2 - start) / range_len
        plot_end = plot_start + 0.1
        mate_chrom = bp1_chrom
        mate_pos = bp1
      }
      pushViewport(viewport(
        x = plot_start,
        y = 1 / 3,
        width = (plot_end - plot_start),
        height = 0.15,
        just = c('left', 'top'),
        clip = 'off'
      ))
      tra_label_pos_y = 4 + runif(1, 0, 4)
      tra_tail_pos_y = -7
      grid.lines(
        x = c(0, 0),
        y = c(tra_label_pos_y, tra_tail_pos_y),
        gp = gpar(col =  sv_color[['TRA']], lwd = 0.2)
      )
      
      if (direction == '5to5'){
       bottom_sign = 1
       top_sign = 1
      }else if(direction == '3to3'){
        bottom_sign = -1 
        top_sign = -1 
      }else if(direction %in% c('3to5', '5to3')){
        if(direction == '3to5'){
          bottom_sign = -1
          top_sign = 1
        }else{
          bottom_sign = 1
          top_sign = -1
        }
        
        if(chrom!=bp1_chrom){
          bottom_sign = -bottom_sign
          top_sign = -top_sign
        }
      }
      
      # print(paste0(top_sign,  ' ', bottom_sign))
      if(top_sign==1){
        text_just = 'left'
        text_placement_scale = 1
      }else{
        text_just = 'right'
        text_placement_scale = 1
      }
      
      grid.lines(
        x = c(0,0.3 * top_sign),
        y = c(tra_label_pos_y, tra_label_pos_y + 0.3),
        gp = gpar(col = sv_color[['TRA']], lwd = 0.1)
        )
      # print(text_just)
      grid.text(paste0(mate_chrom, ':', formatC(mate_pos, format="d", big.mark=",")), x = 0.3*top_sign*text_placement_scale, y = tra_label_pos_y + 0.3, gp=gpar(cex=1, fontsize=6), just=text_just)
      grid.lines(
        x = c(0,0.1 * bottom_sign),
        y = c(tra_tail_pos_y,tra_tail_pos_y - 0.1),
        gp = gpar(col = sv_color[['TRA']], lwd = 0.1)
        )
      upViewport()
    }
  }



sv_curve_vec = Vectorize(sv_curve)

# GET mutation specific copynumber
parse_tumor_cell_fraction<-function(abscn_path){
  abscn_file = file(abscn_path, open='r')
  headers = readLines(abscn_file, n = 3)
  close(abscn_file)
  as.double(str_split(string = headers[3], pattern='\t')[[1]][2])
}
get_regional_abscn<-function(sample, chrom, position, abscn_path){
  abscn_df  = as_tibble(read.table(
    file=abscn_path,
    skip = '#',
    col.names = c('chromosome', 'pos', 'tumordepth', 'normaldepth', 'abscn')
  )) #, delim=' ', , na='NA')
  abscn_region_df = abscn_df %>% filter(chromosome==chrom) %>% filter(between(position-pos, 0, 100000))
  abscn = abscn_region_df$abscn[1]
  abscn
}
get_regional_abscn_vec = Vectorize(get_regional_abscn)
get_mutated_copynumber<- function(raw_vaf, total_abscn, tumor_cell_fraction){
  raw_vaf*(2/tumor_cell_fraction + (total_abscn -2))
}
get_mutated_copynumber_vec = Vectorize(get_mutated_copynumber)

#################################
# PLOT

plot_sv<-function(sample, chrom, start, end, snv_vaf_path, abscn_path, sv_path, output_dir, save=F, format='pdf'){
  # Data to use
  # snv_data = paste0('/Users/cyoon/Dropbox/Research/julab/myeloma/sv_reconstruction/', sample, '_mns.pon.filtered.vep.vcf.gz.vaf.txt')
  # abscn_data = paste0('/Users/cyoon/Dropbox/Research/julab/myeloma/sv_reconstruction/', sample, '.100kbcov.absCN')
  # sv_data = paste0('/Users/cyoon/Dropbox/Research/julab/myeloma/sv_reconstruction/', sample, '_sv_filtered.intersectFilter.vcf.svanno.txt')
  new_folders_to_make = paste0(output_dir, '/', c('kataegis', 'snv_df', 'figure', 'status'))
  for(afolder in new_folders_to_make){
    system(paste0('mkdir -p ', afolder))
  }
  print('plot_sv function starting')

  # Set figure size. 
  # 25mb -> 1 inch width
  
  # overal plot margins (in inches)
  left_outer_margin =  0.25
  right_outer_margin = 0.25
    
  # plotViewport margins (in lines)
  bottom_margin = 2
  left_margin = 3
  top_margin = 0
  right_margin = 1
    
  range_len = end - start
  figure_scaling_factor = 25000000
  figure_width = range_len/figure_scaling_factor + (left_margin + right_margin)/5  + left_outer_margin + right_outer_margin# 5 lines = 1 inch
  
  if(end == chrom_sizes[[chrom]]){
    region_string=''
  }else{
    region_string=paste0(':', start, '-', end)
  } 
  
  if(save==F){
    grid.newpage()
  }else{
    figure_name = paste0(sample,'_chr', chrom, gsub(pattern = ':', replacement = '_', x=region_string), '.', format)
    if(format=='pdf'){
      print(paste0(output_dir, "/figure/", figure_name))
      pdf(paste0(output_dir, "/figure/", figure_name), width = figure_width, height=7, useDingbats = F)
    }else if(format=='png'){
      png(paste0(output_dir, "/figure/", figure_name))
    }
  }
  title_height = 0.15
  n_subplot = 5
  subplot_heights = (1-title_height)/n_subplot
  
  
  # create a viewport to allow some margins
  width_fig = range_len/figure_scaling_factor + (left_margin + right_margin)/5 
  pushViewport(viewport(x=unit(left_outer_margin, 'inch'), y=unit(0, 'npc'), width=unit(width_fig, 'inch'), height=unit(1, 'npc'), just=c('left', 'bottom')))

    # Figure Title
  pushViewport(viewport(
    x = 0,
    y = 1-title_height,
    height = title_height,
    just = c('left', 'bottom')
  ))
  pushViewport(plotViewport(c(bottom_margin, left_margin, top_margin, right_margin)))
  grid.text(paste0(sample, ' chr', chrom, region_string), x=-0.025, y=1, just=c('left','top'))
  
  popViewport(2)
  
  #################################
  # Ideogram
  ideogram_height = 0.02 
  pushViewport(viewport(
    x=0, 
    y=subplot_heights*4, 
    height=ideogram_height, 
    just=c('left', 'bottom')))
  pushViewport(plotViewport(c(0, left_margin, 0, right_margin)))
  draw_ideogram(paste0('chr', chrom), start,  end, end-start, xpos=0, ypos=0.5, height=1, ref='grch37')
  
  popViewport(2)
  #################################
  # ABSCN
  print('abscn plotting')
  pushViewport(viewport(
    x = 0,
    y = subplot_heights * 3,
    height = subplot_heights,
    just = c('left', 'bottom')
  ))
  pushViewport(plotViewport(c(bottom_margin, left_margin, top_margin, right_margin)))
  abscn_df  = as_tibble(read.table(
    abscn_path,
    skip = '#',
    col.names = c('chromosome', 'pos', 'tumordepth', 'normaldepth', 'abscn')
  )) #, delim=' ', , na='NA')
  abscn_roi = abscn_df %>% filter(chromosome == chrom) %>% filter(between(pos, start, end))
  x = (abscn_roi$pos - start) / range_len
  yrange = range(abscn_roi$abscn, na.rm = T)
  max_y = max(4.1, ceiling(yrange[2]))
  y = abscn_roi$abscn / max_y
  
  if(dim(abscn_roi)[1]>=1){
    grid.points(
      x = x,
      y = y,
      pch = 16, size=unit(0.3, 'char')
    )  
  }
  for(i in 1:floor(max_y/2)){
    grid.lines(x=c(0, max(x)), y=c(2*i/max_y, 2*i/max_y), gp=gpar(lty=2, alpha=0.3))
    
  }
  # grid.lines(x=c(0, max(x)), y=c(4/max_y, 4/max_y), gp=gpar(lty=2, alpha=0.3))
  axis_unit = 20000000 # tick mark every 20mb #range_len / 5
  axis_log = nchar(trunc(round(axis_unit))) - 1
  axis_unit_rounded = round(axis_unit / 10 ** (axis_log)) * 10 ** axis_log
  axis_scaling_factor = range_len / axis_unit_rounded
  xaxis_actual_pos = seq(start, end + axis_unit_rounded, by = axis_unit_rounded)
  xaxis_plot_pos = (xaxis_actual_pos - start) / range_len
  xaxis_label_text = simplify_coordinate_vec(xaxis_actual_pos)
  grid.xaxis(at = xaxis_plot_pos, label = xaxis_label_text)
  grid.yaxis(at = 0:max_y / max_y, label = paste0(0:max_y))
  grid.text('absCN', x=unit(-2.5, 'line'), rot=90)
  if(chrom=='22'){
    igl_start = 22385331
    igl_end = 23265206
    grid.rect(x=igl_start/range_len, y=1.1, width=(igl_end-igl_start)/range_len, height=0.1, gp=gpar(fill='blue', alpha=0.4))
    grid.text('IGL', x=igl_start/range_len, y=1.25, gp=gpar(fontface=3))
  }
  
  
  if(chrom=='2'){
    igk_start = 89156673
    igk_end = 90274235
    grid.rect(x=igk_start/range_len, y=1.1, width=(igk_end-igk_start)/range_len, height=0.1, gp=gpar(fill='blue', alpha=0.4))
    grid.text('IGK', x=igk_start/range_len, y=1.25, gp=gpar(fontface=3))
  }
  
  
  if(chrom=='14'){
    igh_start = 106053225
    igh_end = 107283280
    grid.rect(x=igh_start/range_len, y=1.1, width=(igh_end-igh_start)/range_len, height=0.1, gp=gpar(fill='blue', alpha=0.4))
    grid.text('IGH', x=igh_start/range_len, y=1.25, gp=gpar(fontface=3))
  }
  
  popViewport(2)
  #################################
  # SV
  print('plotting SV')
  sv_df = as_tibble(read.table(
    sv_path, header=T))
  
  # split bp into chrom, start, end
  sv_df_long = sv_df %>% separate(bp1, into = c('bp1_chrom', 'bp1_start', 'bp1_end')) %>% separate(bp2, into =
                                                                                                     c('bp2_chrom', 'bp2_start', 'bp2_end'))
  sv_df_long$bp1_start = as.integer(sv_df_long$bp1_start)
  sv_df_long$bp2_start = as.integer(sv_df_long$bp2_start)
  # find SVs that are in this region
  sv_df_roi = sv_df_long %>% filter((bp1_chrom == chrom &
                                       between(bp1_start, start, end)) |
                                      (bp2_chrom == chrom & between(bp1_start, start, end)))
  sv_df_roi$svtype = as.character(sv_df_roi$svtype)
  sv_df_roi$direction = as.character(sv_df_roi$direction)
  pushViewport(viewport(
    x = 0,
    y = subplot_heights * 4 + ideogram_height,
    height = subplot_heights,
    just = c('left', 'bottom'), 
    clip = 'off'
  ))
  pushViewport(plotViewport(c(0, left_margin, 0, right_margin), clip = F))
  grid.lines(x = c(0, 1),
             y = 1 / 3,
             gp = gpar(col = 'gray'))
  grid.lines(x = unit(c(-0.3, 0), 'inches'),
             y = 1 / 3,
             gp = gpar(col = 'gray'))
  grid.text(
    'DEL',
    x = unit(-0.3, 'inches'),
    y = unit(2 / 3, 'native') + unit(1, 'mm'),
    just = c('left', 'bottom'),
    gp = gpar(col = 'gray')
  )
  grid.text(
    'DUP',
    x = unit(-0.3, 'inches'),
    y = unit(2 / 3, 'native') - unit(1, 'mm'),
    just = c('left', 'top'),
    gp = gpar(col = 'gray')
  )
  grid.lines(x = c(0, 1),
             y = 2 / 3,
             gp = gpar(col = 'gray'))
  grid.lines(x = unit(c(-0.3, 0), 'inches'),
             y = 2 / 3,
             gp = gpar(col = 'gray'))
  grid.text(
    'HH',
    x = unit(-0.3, 'inches'),
    y = unit(1 / 3, 'native') + unit(1, 'mm'),
    just = c('left', 'bottom'),
    gp = gpar(col = 'gray')
  )
  grid.text(
    'TT',
    x =unit(-0.3, 'inches'),
    y = unit(1 / 3, 'native') - unit(1, 'mm'),
    just = c('left', 'top'),
    gp = gpar(col = 'gray')
  )
  if(dim(sv_df_roi)[1]>=1){
    for (i in 1:dim(sv_df_roi)[1]) {
      # print(i)
      sv_curve(
        sv_df_roi$bp1_chrom[i],
        sv_df_roi$bp1_start[i],
        sv_df_roi$bp2_chrom[i],
        sv_df_roi$bp2_start[i],
        start,
        range_len,
        sv_df_roi$svtype[i],
        sv_df_roi$direction[i]
      )
    }
    
    
    # print(as.character(sv_df_roi[i, 'svtype']))
    # cat(as.integer(sv_df_roi[i, 'bp1_start']), as.integer(sv_df_roi[i, 'bp2_start']), start, range_len, as.character(sv_df_roi[i, 'svtype']), as.character(sv_df_roi[i, 'direction']))
  }
  upViewport(2)
  #################################
  # SNV
  print('plotting SNVs')
  snv_df = as_tibble(read.table(snv_vaf_path, header = T))
  # print(snv_df)
  snv_df_roi = snv_df %>% filter(chromosome==chrom) %>% filter(pos < end) %>% filter(pos > start)
  substitution_color= c('C>G'= 'black', 'C>A'= 'blue', 'C>T'= 'red', 'T>A'= 'magenta', 'T>C'= 'yellow', 'T>G'= 'green')
  # substitution_color= c('C>G'= 'white', 'C>A'= 'white', 'C>T'= 'red', 'T>A'= 'white', 'T>C'= 'white', 'T>G'= 'white')
  # substitution_color= c('C>G'= NULL, 'C>A'= NULL, 'C>T'= 'red', 'T>A'= NULL, 'T>C'= NULL, 'T>G'= NULL)
  # print(snv_df_roi)
  if(dim(snv_df_roi)[1]>0){
    snv_df_roi$sample = sample
    snv_df_roi_abscn = snv_df_roi %>% mutate(abscn=get_regional_abscn_vec(sample, chrom, pos, abscn_path))
    print(abscn_path)
    print(parse_tumor_cell_fraction(abscn_path))
    snv_df_roi_abscn$tumor_cell_fraction = parse_tumor_cell_fraction(abscn_path)
    snv_df_roi_abscn_mutcopy = snv_df_roi_abscn %>% mutate(mutcopy = get_mutated_copynumber_vec(vaf, abscn, tumor_cell_fraction))
    # write_delim(snv_df_roi_abscn_mutcopy, paste0('~/Dropbox/Research/julab/myeloma/figures/snv_df/', sample, '_', chrom, '_snv_df.txt'), delim = '\t')
    substitution_x = (snv_df_roi$pos - start)/range_len
  }
  
  pushViewport(viewport(
    x = 0,
    y = subplot_heights*2,
    height = subplot_heights,
    just = c('left', 'bottom')
  ))
  # grid.rect() 
  pushViewport(plotViewport(c(bottom_margin, left_margin, top_margin, right_margin)))
  grid.xaxis(at = xaxis_plot_pos, label = xaxis_label_text)
  grid.yaxis(gp = gpar(fontsize = 8))
  grid.text('VAF', x=unit(-2.5, 'line'), rot=90)
  if(dim(snv_df_roi)[1]>=1){
    grid.points(x=substitution_x, y=snv_df_roi$vaf, pch=16, size=unit(1, 'mm'), gp=gpar(col=substitution_color[snv_df_roi$substitution]))
  }
  
  ######### Mutated Copy
  popViewport(2)
  pushViewport(viewport(
    x = 0,
    y = subplot_heights*1,
    height = subplot_heights,
    just = c('left', 'bottom')
  ))
  pushViewport(plotViewport(c(bottom_margin, left_margin, top_margin, right_margin)))
  
  grid.xaxis(at = xaxis_plot_pos, label = xaxis_label_text)
  
  
  if(dim(snv_df_roi)[1]>=1){
    grid.points(x=substitution_x, y=snv_df_roi_abscn_mutcopy$mutcopy/max_y, pch=16, size=unit(1, 'mm'), gp=gpar(col=substitution_color[snv_df_roi$substitution]))
    grid.yaxis(at=0:max_y/max_y, label=0:max_y, gp = gpar(fontsize = 8))
    for(i in 1:floor(max_y/2)){
      grid.lines(x=c(0, max(x)), y=c(2*i/max_y, 2*i/max_y), gp=gpar(lty=2, alpha=0.3))
      
    }  
  }else{
    # place holder for the axis when no data points are available
    max_y=2
    grid.yaxis(at=0:max_y/max_y, label=0:max_y, gp = gpar(fontsize = 8))
  }
  
  grid.text('Mutated CN', x=unit(-2.5, 'line'), rot=90)
  
  #################################
  # Mutation Distance
  print('plotting rainfall')
  popViewport(2)
  pushViewport(viewport(
    x = 0,
    y = subplot_heights*0,
    height = subplot_heights,
    just = c('left', 'bottom')
  ))
  pushViewport(plotViewport(c(bottom_margin, left_margin, top_margin, right_margin)))
  if(dim(snv_df_roi)[1] > 0){
    mut_distance = snv_df_roi_abscn_mutcopy %>% mutate(distance=log10(pos-lag(pos, default = 0)))
    mut_distance$index = seq(1:dim(mut_distance)[1])
    mut_distance$kataegis = NA
    mut_distance_under1000 = mut_distance %>% filter(distance<=3)
    k = split(mut_distance_under1000$index, cumsum(c(1, diff(mut_distance_under1000$index) != 1)))
    kataegis_count = 0
    
    # Definition of kataegis: 
    for(i in names(k)){
      if(length(k[[i]])>=6){
        kataegis_count = kataegis_count + 1
        for(index in k[[i]]){
          mut_distance$kataegis[index] = paste0(sample, '_', chrom,'_k', as.character(kataegis_count))
        }
      }
    }
    kataegis_foci = mut_distance %>% filter(!is.na(kataegis)) %>% group_by(kataegis) %>% summarise(start=min(pos), end=max(pos))
    write_delim(mut_distance, paste0(output_dir,  '/snv_df/', sample, '_', chrom, '_snv_df.txt'), delim = '\t')
    grid.points(x=substitution_x, y=mut_distance$distance/8, pch=16, size=unit(1, 'mm'), gp=gpar(col=substitution_color[mut_distance$substitution]))
    grid.yaxis(at=0:8/8, label=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression(10^7), expression(10^8)), gp = gpar(fontsize = 8))
    
    if(kataegis_count >=1){
      kataegis_snvs = mut_distance%>% filter(!is.na(kataegis))
      kataegis_snvs$substitution_x = (kataegis_snvs$pos - start)/range_len
      write_delim(kataegis_foci, paste0(output_dir, '/kataegis/', sample, '_', chrom, '_kataegis.txt'), delim = '\t')
      grid.points(x=kataegis_snvs$substitution_x, y=rep(-0.15, dim(kataegis_snvs)[1]), pch=24, gp=gpar(fill='yellow'), size=unit(0.03, 'npc'))
    }
  }
  else{
    # place holder for the axis when no data points are available
    max_y=2
    grid.yaxis(at=0:8/8, label=c(expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression(10^7), expression(10^8)), gp = gpar(fontsize = 8))
  }
  
  
  grid.xaxis(at = xaxis_plot_pos, label = xaxis_label_text)
  
  grid.text('Mutation distance (bp)', x=unit(-4.5, 'line'), rot=90, gp=gpar(fontsize=7))
  
  
  if(save==T){
    dev.off()
  }
  
} 



############################
if (!interactive()) {
  
  # Plot all chromosomes for all samples in MM data
  for(chrom in names(chrom_sizes)){
    print(paste0('chrom: ', chrom))
    print(paste0('chromsize: ', chrom_sizes[[chrom]]))
    print(paste0('snv_vaf_path: ', snv_vaf_path))
    print(paste0('abscn_path: ', abscn_path))
    print(paste0('svpath: ', sv_path))
    print(paste0('outputdir: ', output_dir))
    print(paste0('sample_name: ', sample_name))
    # print(sample_name, chrom, 0, chrom_sizes[[chrom]], snv_vaf_path, abscn_path, sv_path, output_dir)
    # (sample, chrom, start, end, snv_vaf_path, abscn_path, sv_path, output_dir, save=F, format='pdf')
    plot_sv(sample_name, chrom, 0, chrom_sizes[[chrom]], snv_vaf_path, abscn_path, sv_path, output_dir, save=T, format='pdf')
  }
  
  # write empty file to mark job completion
  x <- data.frame()
  write.table(x, file=paste0(output_dir, '/status/', sample_name, '.svplot.done'), col.names=FALSE)
  
  
}
