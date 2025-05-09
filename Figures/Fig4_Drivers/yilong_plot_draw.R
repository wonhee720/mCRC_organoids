# args[1] = path/to/100kbcov.absCN
# args[2] = path/to/sv.list
# args[3] = chromosome of interest
# args[4] = (optional) start position
# args[5] = (optional) end position - required if supplying start position
# args[6] = (optional) SV cluster number
# args[7] = (optional) branch among the cluster (SVs accumulated to that branch from the root)
# Modify the source function directory in order to get mosdepth or sequenza as input - LWH 21.12.23 

suppressMessages(library(tidyverse))

#source("/home/users/wonhee720/tools/scripts/SV/yilong_plot/yilong_plot_function_mosdepth_gc.R")
source("/home/users/wonhee720/tools/scripts/SV/yilong_plot/yilong_plot_function_mosdepth_gc.P57.R")

# Use Patient specific .R yilong plot functions are ready (ex-"/home/users/wonhee720/tools/scripts/SV/yilong_plot/yilong_plot_function_mosdepth_gc.69-LS78-O5.R", "/home/users/wonhee720/tools/scripts/SV/yilong_plot/yilong_plot_function_mosdepth_gc.P57.R")


args <- commandArgs(trailingOnly = T)

print(length(args))

if (length(args) == 3){
  draw_svsketch(args[1], args[2], args[3], save_plot = T)
  } else if (length(args) == 6) {     # <- receives cluster information
  draw_svsketch(args[1], args[2], args[3], args[4], args[5], args[6], save_plot = T)
  } else if (length(args) == 7) {     # <- receives cluster & branch information
  draw_svsketch(args[1], args[2], args[3], args[4], args[5], args[6], args[7], save_plot = T)
  } else {
  draw_svsketch(args[1], args[2], args[3], args[4], args[5], save_plot = T)
}
