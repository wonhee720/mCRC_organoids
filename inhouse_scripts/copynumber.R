library(tidyverse)
library(cowplot)
library(mclust)
suppressMessages(library(grid))
suppressMessages(library(optparse))
suppressMessages(library(grideogram))
# gridegram is not available 


setwd("/home/users/kimin/projects/03_Mucosal_Melanoma/")
list.files()

reference = "/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta"
abscn_path = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/MM001-Tumor-Blood.100kbcov.absCN"
output_dir = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/copynumber/"

#### General function setup 
parse_chrom_size <- function(reference_path){
  fasta_index = paste0(reference_path, '.fai')
  reference=read_delim(fasta_index, delim='\t', col_names=c('chromosome', 'size', 'a', 'b', 'c'))
  structure(reference$size, names=reference$chromosome)
}

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

chrom_sizes = parse_chrom_size(reference)
chrom_sizes = chrom_sizes[str_detect(names(chrom_sizes), pattern = "^[0-9XY]*$")] # only get autosomal + X Y chromosomes, ignore small contigs



#### GET mutation specific copynumber
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



#### Plot function
plot_cn<-function(sample, chrom, start, end, abscn_path, output_dir, format='pdf'){
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
  
  # Generate variable called region_string
  if(end == chrom_sizes[[chrom]]){
    region_string=''
  }else{
    region_string=paste0(':', start, '-', end)
  }
  
  # File name adjusting
  figure_name = paste0(sample,'_chr', chrom, gsub(pattern = ':', replacement = '_', x=region_string), '.', format)
  if(format=='pdf'){
    print(paste0(output_dir, figure_name))
    pdf(paste0(output_dir, figure_name), width = figure_width, height=7, useDingbats = F)
  }else if(format=='png'){
    png(paste0(output_dir, figure_name))
  }
  
  title_height = 0.15
  n_subplot = 1
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
  
  #### ABSCN
  pushViewport(viewport(
    x = 0,
    y = subplot_heights * 1,
    height = subplot_heights,
    just = c('left', 'bottom')
  ))
  pushViewport(plotViewport(c(bottom_margin, left_margin, top_margin, right_margin)))
  abscn_df  = as_tibble(read.table(
    abscn_path,
    skip = '#',
    col.names = c('chromosome', 'pos', 'tumordepth', 'normaldepth', 'abscn')
  )) #, delim=' ', , na='NA'
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
  #grid.lines(x=c(0, max(x)), y=c(4/max_y, 4/max_y), gp=gpar(lty=2, alpha=0.3))
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
  
  dev.off()
  
}
  
  
  



#### RUN SCRIPT

abscn_path = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/MM001-Tumor-Blood.100kbcov.absCN"
abscn_path = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/MM001-Tumor-Normal.100kbcov.absCN"
abscn_path = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/MM001-Normal-Blood.100kbcov.absCN"


sample = 'Tumor-Blood'
start = 0
end = chrom_sizes[[1]]
chrom = 1

for(chrom in names(chrom_sizes)){
  plot_cn(sample, 1, 0, chrom_sizes[[1]], abscn_path, output_dir, format='pdf')
}





#### Test 
sample = "Tumor"
average_cov = 70.29

single_cov = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/MM001-Tumor.mpileup.100kbcov"
single_cov_df  = as_tibble(read.table(
  single_cov,
  skip = '#',
  header = TRUE,
  col.names = c('chromosome', 'pos', 'N', 'effective', 'GC', 'sum', 'avg')
))

pdf(paste0(output_dir, figure_name), width=14, height=7, useDingbats = F)

single_cov_roi = single_cov_df %>% filter(chromosome == 1) %>% filter(between(pos, 0, chrom_sizes[[1]]))
range_len = chrom_sizes[[1]]
x = (single_cov_roi$pos) / range_len
y = single_cov_roi$avg
y[y >= average_cov * 4] = average_cov*4
yrange = range(y, na.rm = T)
max_y = max(4.1, ceiling(yrange[2]))
y = single_cov_roi$avg / max_y

# overal plot margins (in inches)
left_outer_margin =  0.25
right_outer_margin = 0.25

# plotViewport margins (in lines)
bottom_margin = 0.1
left_margin = 0.1
top_margin = 0.1
right_margin = 0.1

title_height = 0.05

width_fig = range_len/figure_scaling_factor + 
  (left_margin + right_margin)/5 

figure_width = range_len/figure_scaling_factor + 
  (left_margin + right_margin)/5  + 
  left_outer_margin + right_outer_margin # 5 lines = 1 inch

pushViewport(viewport(
  x=0, 
  y=0, 
  width=1, 
  height=1, 
  just=c('left', 'bottom')))

grid.rect(gp=gpar(lty="dashed"))

pushViewport(viewport(
  x = 0,
  y = 1-title_height,
  height = title_height,
  just = c('left', 'bottom')
))

grid.rect(gp=gpar(lty="dashed"))

#pushViewport(plotViewport(
#  c(bottom_margin, left_margin, top_margin, right_margin)))

#grid.rect(gp=gpar(lty="dashed"))

grid.text(paste0(sample, ' Chr ', chrom), x=0.1, y=0.5, just=c('left','center'))

popViewport(2)

pushViewport(viewport(
  x = 0,
  y = 0,
  height = 1 - title_height,
  just = c('left', 'bottom')
))

pushViewport(viewport(
  x = 0.5,
  y = 0.5,
  height = 0.8,
  width = 0.8
  ))

grid.rect(gp=gpar(lty="dashed"))


#pushViewport(plotViewport(
#  c(bottom_margin, left_margin, top_margin, right_margin)))

grid.points(
  x = x,
  y = y,
  pch = 16, size=unit(0.3, 'char')
  )

grid.lines(x=c(0,max(x)), y=average_cov/max_y, gp=gpar(col = 'red', lwd=2))

axis_unit = 20000000 
axis_scaling_factor = range_len / axis_unit
xaxis_actual_pos = c(seq(0, end, by = axis_unit), end)
xaxis_plot_pos = (xaxis_actual_pos) / range_len
xaxis_label_text = simplify_coordinate_vec(xaxis_actual_pos)
grid.xaxis(at = xaxis_plot_pos, label = xaxis_label_text)
yaxis_plot_pos = seq(0, 1+average_cov/max_y, by = average_cov/max_y)
grid.yaxis(at = yaxis_plot_pos, label = paste0(0:4))
grid.text('copy number from coverage', x=-0.1, rot=90)

 popViewport(2)

 dev.off()

 