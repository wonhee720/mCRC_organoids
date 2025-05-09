suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
#suppressMessages(library(mclust))
suppressMessages(library(grid))
suppressMessages(library(optparse))
#suppressMessages(library(grideogram))

## sample table used for iteration
Blood = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/depth/MM001-Blood.mpileup.100kbcov"
Normal = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/depth/MM001-Normal.mpileup.100kbcov"
Tumor = "/home/users/kimin/projects/03_Mucosal_Melanoma/02_vcf/depth/MM001-Tumor.mpileup.100kbcov"
Blood_average_cov = 38.88
Normal_average_cov = 58.76
Tumor_average_cov = 70.52
sample = rbind(c('Blood', 'Normal', 'Tumor'), c(Blood, Normal, Tumor), c(Blood_average_cov, Normal_average_cov, Tumor_average_cov))

## initial variables and functions to use later
setwd("/home/users/kimin/projects/03_Mucosal_Melanoma/")
output_dir = "/home/users/kimin/projects/03_Mucosal_Melanoma/06_fig/copynumber_post_fastp/"
reference = "/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta"


### GC Check
sample_path = Blood
sample_data = as_tibble(read.table(
  sample_path,
  skip = '#',
  header = TRUE,
  col.names = c('chromosome', 'pos', 'N', 'effective', 'GC', 'sum', 'avg')
))

sample_data = sample_data %>% filter(chromosome %in% names(chrom_sizes))
sample_data$pos_genome = sample_data$pos+chrom_sizes_acc[as.character(sample_data$chromosome)]

ggplot(data = sample_data, aes(x=avg, y=GC))+geom_point()+xlim(0,200)+ylim(0,1)


xrange = max(sample_data$pos_genome)
x = sample_data$pos_genome / xrange
y = sample_data$GC
yrange = range(y, na.rm = T)
y = y / yrange

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

grid.text(paste0(sample_name, ' Whole Genome'), x=(0.2/6), y=0.5, just=c('left','center'))

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
  width = 1 - (0.2/3) 
))

grid.rect(gp=gpar(lty="dashed"))


xaxis_actual_pos = c(0)
for(i in c(1:24)){xaxis_actual_pos[i + 1] = sum(chrom_sizes[1:i])}
xaxis_plot_pos = xaxis_actual_pos / xrange
grid.xaxis(at = xaxis_plot_pos, label = FALSE)
yaxis_plot_pos = seq(0, 1, by = (sample_average_cov/2)/yrange)
grid.yaxis(at = yaxis_plot_pos, label = paste0(0:(length(yaxis_plot_pos)-1)))
for(i in c(1:24)){
  grid.text(names(chrom_sizes)[i], x= (xaxis_plot_pos[i+1] + xaxis_plot_pos[i])/2, y=-0.05)
}
grid.text('copy number from coverage', x=-0.03, rot=90)
grid.text('chromosome', y=-0.1)

for(i in c(1:24)){
  if(i %% 2 == 1){bgcolor = 8}
  else{bgcolor = FALSE}
  grid.rect(x = xaxis_plot_pos[i], y = 0, just = c('left','bottom'), height = 1, width = xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE, alpha = 0.2))
}

grid.points(
  x = x,
  y = y,
  pch = 16, size=unit(0.3, 'char')
)

grid.lines(x=c(0,max(x)), y=sample_average_cov/yrange, gp=gpar(col = 'red', lwd=2))

popViewport(2)





bottom_margin = 0.1
left_margin = 0.1
top_margin = 0.1
right_margin = 0.1

title_height = 0.05

parse_chrom_size <- function(reference_path){
  fasta_index = paste0(reference_path, '.fai')
  reference=read_delim(fasta_index, delim='\t', col_names=c('chromosome', 'size', 'a', 'b', 'c'))
  structure(reference$size, names=reference$chromosome)
}
chrom_sizes = parse_chrom_size(reference)
chrom_sizes = chrom_sizes[str_detect(names(chrom_sizes), pattern = "^[0-9XY]*$")] # only get autosomal + X Y chromosomes, ignore small contigs
chrom_sizes_acc= c()
chrom_sizes_acc[1] = 0
for(i in c(2:24)){
  chrom_sizes_acc[i] = sum(chrom_sizes[1:i - 1])
}
chrom_sizes_acc = structure(chrom_sizes_acc, names = names(chrom_sizes))

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



## Plot function
plot<-function(sample, chrom, start, end, data, average_cov, output_dir){
  figure_name = paste0(sample,'_chr_', chrom, '.pdf')
  pdf(paste0(output_dir, figure_name), width=14, height=7, useDingbats = F)
  
  data_chrom = data %>% filter(chromosome == chrom) %>% filter(between(pos, start, end))
  range_len = end - start
  x = data_chrom$pos / range_len
  y = data_chrom$avg
  y[y >= average_cov * 4] = average_cov*4
  yrange = range(y, na.rm = T)
  yrange = max(average_cov*4, ceiling(yrange[2]))
  y = y / yrange
  
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
  
  grid.text(paste0(sample, ' Chr ', chrom), x=0.05, y=0.5, just=c('left','center'))
  
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
    width = 0.9
  ))
  
  grid.rect(gp=gpar(lty="dashed"))
  
  grid.points(
    x = x,
    y = y,
    pch = 16, size=unit(0.3, 'char')
  )
  
  grid.lines(x=c(0,max(x)), y=average_cov/yrange, gp=gpar(col = 'red', lwd=2))
  
  axis_unit = 20000000 
  axis_scaling_factor = range_len / axis_unit
  xaxis_actual_pos = c(seq(0, chrom_sizes[[chrom]], by = axis_unit), chrom_sizes[[chrom]])
  xaxis_plot_pos = (xaxis_actual_pos) / range_len
  xaxis_label_text = simplify_coordinate_vec(xaxis_actual_pos)
  grid.xaxis(at = xaxis_plot_pos, label = xaxis_label_text)
  yaxis_plot_pos = seq(0, 1, by = (average_cov/2)/yrange)
  grid.yaxis(at = yaxis_plot_pos, label = paste0(0:(length(yaxis_plot_pos)-1)))
  grid.text('copy number from coverage', x=-0.1, rot=90)
  
  popViewport(2)
  
  dev.off()
  
}



## Main script
for( i in c(1:3)){
  sample_name = sample[1,i]
  sample_path = sample[2,i]
  sample_average_cov = as.numeric(sample[3,i])
  sample_data = as_tibble(read.table(
    sample_path,
    skip = '#',
    header = TRUE,
    col.names = c('chromosome', 'pos', 'N', 'effective', 'GC', 'sum', 'avg')
  ))
  
  output_dir_name = paste0(output_dir, sample_name, '/')
  system(paste0('mkdir -p ', output_dir_name))
  
  for( chrom in names(chrom_sizes)){
    plot(sample_name, chrom, 0, chrom_sizes[[chrom]], sample_data, sample_average_cov, output_dir_name)
  }
  
}



## Main script 2
for( i in c(1:3)){
  sample_name = sample[1,i]
  sample_path = sample[2,i]
  sample_average_cov = as.numeric(sample[3,i])
  output_dir_name = paste0(output_dir, sample_name, '/')
  figure_name = paste0(sample_name,'_overview', '.pdf')
  system(paste0('mkdir -p ', output_dir_name))
  pdf(paste0(output_dir_name, figure_name), width=21, height=7, useDingbats = F)
  
  sample_data = as_tibble(read.table(
    sample_path,
    skip = '#',
    header = TRUE,
    col.names = c('chromosome', 'pos', 'N', 'effective', 'GC', 'sum', 'avg')
  ))
  
  sample_data = sample_data %>% filter(chromosome %in% names(chrom_sizes))
  sample_data$pos_genome = sample_data$pos+chrom_sizes_acc[as.character(sample_data$chromosome)]
  
  xrange = max(sample_data$pos_genome)
  x = sample_data$pos_genome / xrange
  y = sample_data$avg
  y[y >= sample_average_cov * 4] = sample_average_cov*4
  yrange = range(y, na.rm = T)
  yrange = max(sample_average_cov*4, ceiling(yrange[2]))
  y = y / yrange
  
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
  
  grid.text(paste0(sample_name, ' Whole Genome'), x=(0.2/6), y=0.5, just=c('left','center'))
  
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
    width = 1 - (0.2/3) 
  ))
  
  grid.rect(gp=gpar(lty="dashed"))
  
  
  xaxis_actual_pos = c(0)
  for(i in c(1:24)){xaxis_actual_pos[i + 1] = sum(chrom_sizes[1:i])}
  xaxis_plot_pos = xaxis_actual_pos / xrange
  grid.xaxis(at = xaxis_plot_pos, label = FALSE)
  yaxis_plot_pos = seq(0, 1, by = (sample_average_cov/2)/yrange)
  grid.yaxis(at = yaxis_plot_pos, label = paste0(0:(length(yaxis_plot_pos)-1)))
  for(i in c(1:24)){
    grid.text(names(chrom_sizes)[i], x= (xaxis_plot_pos[i+1] + xaxis_plot_pos[i])/2, y=-0.05)
  }
  grid.text('copy number from coverage', x=-0.03, rot=90)
  grid.text('chromosome', y=-0.1)
  
  for(i in c(1:24)){
    if(i %% 2 == 1){bgcolor = 8}
    else{bgcolor = FALSE}
    grid.rect(x = xaxis_plot_pos[i], y = 0, just = c('left','bottom'), height = 1, width = xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE, alpha = 0.2))
  }
  
  grid.points(
    x = x,
    y = y,
    pch = 16, size=unit(0.3, 'char')
  )
  
  grid.lines(x=c(0,max(x)), y=sample_average_cov/yrange, gp=gpar(col = 'red', lwd=2))
  
  popViewport(2)
  
  dev.off()
}
