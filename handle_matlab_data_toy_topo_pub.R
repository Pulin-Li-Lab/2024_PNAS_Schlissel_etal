setwd('/path/to/working/directory') ## this is where the plots will be save

library(readr)
library(dplyr)
library(tidyr)
library(R.matlab)
library(ggplot2)
library(wesanderson)


build_msd_data = function(data , old_data = NULL , cell_size = 0 , fracFree = 0.05 , fastRate = 20 , slowRate = 0.5 , topo = 'square'){
  new_dat = NULL
  n_sets = dim(data[[1]])[2] 
  for(i in 1:n_sets){
    this_set = (data[[1]][1,i][[1]][[1]])
    
    to_add = data.frame(time = 1:(unlist(this_set['max.time',1,1])) , 
                        msds = unlist(this_set['msds',1,1]),
                        fluxxed = unlist(this_set['fluxxed',1,1]),
                        model = strsplit(unlist(this_set['name',1,1]),'_')[[1]][1],
                        toTwo = this_set[['matrix',1,1]][1,2] ,
                        toOne = this_set[['matrix',1,1]][2,1] ,
                        cell_size = cell_size,
                        fracFree = fracFree,
                        fastRate = fastRate,
                        slowRate = slowRate , 
                        topo = topo)
    new_dat = rbind(new_dat , to_add)
  }
  
  dat = rbind(old_data , new_dat)
  
  return(dat)
}

read_in_filelist = function(dir , basename , fracFree , cell_size , slowRate = 0.5 , fastRate = 1 , topo = 'square'){
  d = NULL
  for(j in 1:length(basename)){
    f_list = list.files(dir , full.names = T , pattern = basename[j])
    for(i in 1:length(f_list)){
      print(i)
      try({d = rbind(d , build_msd_data(data = readMat(f_list[i])  , fracFree = fracFree , cell_size = cell_size , fastRate = fastRate , topo = topo))})
    }
  }
  return(d)
}

dat = NULL
dir = '/path/to/matlab/output/noOffset'
basename = c('estDRate_toy_agents=5000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=1_condition=constrained')
dat = read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'square')

dir = '/path/to/matlab/output/noOffset'
basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=1_condition=free')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'square'))

dir = '/path/to/matlab/output/squareOffset'
basename = c('estDRate_toy_agents=5000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=1_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'squareOffset'))

dir = '/path/to/matlab/output/squareOffset'
basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=1_condition=free')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'squareOffset'))


## use the range 1000:20000 to calculate the average slope of MSD vs T, i.e. the D_eff
## note the timestep in the simulations is 0.006s, hence the constant below...
range = 1000:20000
dat %>%
  group_by(toOne , toTwo , fastRate, model, fracFree , cell_size , topo) %>%
  summarize(d_eff = lm(msds[range] ~ time[range])$coefficients[2] / 0.006 / 4) %>% ## MSD ~ 4D∆t in 2 dimensions
  group_by(fastRate, model, fracFree , cell_size , topo) %>%
  mutate(d_eff = d_eff , 
         scaled_d_eff = d_eff / d_eff[1]) -> 
  d_eff_dat

unscaled_d = ggplot(d_eff_dat , aes(x = toOne , y = d_eff , col = interaction(topo , model) , group = interaction(topo , model) )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  #coord_cartesian(xlim = c(0,0.5)) +
  ylab('d_eff (µm2/s)') +
  scale_color_manual(values = wes_palette("Zissou1", n = 5 , type='continuous')[c(1,4,2,5)])
unscaled_d

scaled_d = ggplot(d_eff_dat , aes(x = toOne , y = scaled_d_eff , col = interaction(topo , model) , group = interaction(topo , model))) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  #geom_smooth(span = 0.5)+
  #coord_cartesian(xlim = c(0,0.5)) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  ylab('d_eff rel 0%') +
  scale_color_manual(values = wes_palette("Zissou1", n = 5 , type='continuous')[c(1,4,2,5)])
scaled_d

d_eff_plots = gridExtra::arrangeGrob(unscaled_d , scaled_d  , nrow = 1)
#ggsave('./offsetSquare_comp_Deff.pdf' , d_eff_plots , width = length(d_eff_plots) * 5 , height = 3)


unscaled_d = ggplot(d_eff_dat[d_eff_dat$topo == 'square',] , aes(x = toOne , y = d_eff , col = interaction(topo , model) , group = interaction(topo , model) )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  #coord_cartesian(xlim = c(0,0.5)) +
  ylab('d_eff (µm2/s)')+
  scale_color_manual(values = wes_palette("Zissou1", n = 2 , type='continuous'))
unscaled_d

scaled_d = ggplot(d_eff_dat[d_eff_dat$topo == 'square',]  , aes(x = toOne , y = scaled_d_eff , col = interaction(topo , model) , group = interaction(topo , model))) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  #geom_smooth(span = 0.5)+
  #coord_cartesian(xlim = c(0,0.5)) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  ylab('d_eff rel 0%')+
  scale_color_manual(values = wes_palette("Zissou1", n = 2 , type='continuous'))
scaled_d

d_eff_plots = gridExtra::arrangeGrob(unscaled_d , scaled_d  , nrow = 1)
#ggsave('./d_eff_comparison.pdf' , d_eff_plots , width = length(d_eff_plots) * 4.5 , height = 3)


## fastRate comparision figure
dat = NULL
dir = '/path/to/matlab/output/'
basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=0.1_condition=constrained')
dat = read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 0.1 , topo = 'square')

basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=0.5_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 0.5 , topo = 'square'))

basename = c('estDRate_toy_agents=5000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=1_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'square'))

basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=2_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 2 , topo = 'square'))

basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=4_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 4 , topo = 'square'))

basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=10_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 10 , topo = 'square'))


## use the range 1000:20000 to calculate the average slope of MSD vs T, i.e. the D_eff
## note the timestep in the simulations is 0.006s, hence the constant below...range = 1000:20000
dat %>%
  group_by(toOne, fastRate, model, fracFree , cell_size , topo) %>%
  summarize(d_eff = lm(msds[range] ~ time[range])$coefficients[2] / 0.006 / 4) %>% ## MSD ~ 4D∆t in 2 dimensions
  group_by(fastRate, model, fracFree , cell_size , topo) %>%
  mutate(d_eff = d_eff , 
         scaled_d_eff = d_eff / d_eff[1]) -> 
  d_eff_dat

unscaled_d = ggplot(d_eff_dat , aes(x = toOne , y = d_eff , col = log(fastRate) , group = fastRate )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  coord_cartesian(xlim = c(0,1)) +
  ylab('d_eff (µm2/s)') 
unscaled_d

scaled_d = ggplot(d_eff_dat , aes(x = toOne , y = scaled_d_eff , col = log(fastRate) , group = fastRate )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  coord_cartesian(xlim = c(0,1)) +
  ylab('d_eff (rel. tx=0)')
scaled_d

d_eff_plots = gridExtra::arrangeGrob(unscaled_d , scaled_d , nrow = 1)
plot(d_eff_plots)
#ggsave('./d_eff_comparison_by_fastRate.pdf' , d_eff_plots , width = length(d_eff_plots) * 4.5 , height = 3)



## cellSize comparision figure
dat = NULL
dir = '/path/to/matlab/output/'
basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=9_fracFast=0.05_slowRate=0.5_fastRate=1_condition=constrained')
dat = read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 9 , slowRate = 0.5 , fastRate = 1 , topo = 'square')

basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=16_fracFast=0.05_slowRate=0.5_fastRate=1_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 16 , slowRate = 0.5 , fastRate = 1 , topo = 'square'))

basename = c('estDRate_toy_agents=5000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=1_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'square'))

basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=100_fracFast=0.05_slowRate=0.5_fastRate=1_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 100 , slowRate = 0.5 , fastRate = 1 , topo = 'square'))

## use the range 1000:20000 to calculate the average slope of MSD vs T, i.e. the D_eff
## note the timestep in the simulations is 0.006s, hence the constant below...
range = 1000:20000
dat %>%
  group_by(toOne, fastRate, model, fracFree , cell_size , topo) %>%
  #summarize(d_eff = mean(diff(msds[range] , 100) / 100) / 0.006) %>%
  summarize(d_eff = lm(msds[range] ~ time[range])$coefficients[2] / 0.006 / 4) %>% ## MSD ~ 4D∆t in 2 dimensions
  group_by(fastRate, model, fracFree , cell_size , topo) %>%
  mutate(d_eff = d_eff , 
         scaled_d_eff = d_eff / d_eff[1]) -> 
  d_eff_dat

unscaled_d = ggplot(d_eff_dat , aes(x = toOne , y = d_eff , col = factor(cell_size) , group = cell_size )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  coord_cartesian(xlim = c(0,1)) +
  ylab('d_eff (µm2/s)') + 
  scale_color_manual(values = wes_palette("Zissou1", n = 6 , type='continuous'))
unscaled_d

scaled_d = ggplot(d_eff_dat , aes(x = toOne , y = scaled_d_eff , col = factor(cell_size) , group = cell_size )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  coord_cartesian(xlim = c(0,1)) +
  ylab('d_eff (rel. tx=0)') + 
  scale_color_manual(values = wes_palette("Zissou1", n = 6 , type='continuous'))
scaled_d

d_eff_plots = gridExtra::arrangeGrob(unscaled_d , nrow = 1)
#ggsave('./model/d_eff_by_cellSize.pdf' , d_eff_plots , width = length(d_eff_plots) * 4.5 , height = 3)



## fracFree comparision figure
dat = NULL
dir = '/path/to/matlab/output/'
basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.01_slowRate=0.5_fastRate=1_condition=constrained')
dat = read_in_filelist(dir , basename , fracFree = 0.01 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'square')

basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.02_slowRate=0.5_fastRate=1_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.02 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'square'))

basename = c('estDRate_toy_agents=5000_max_time=20000_cellArea=25_fracFast=0.05_slowRate=0.5_fastRate=1_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.05 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'square'))

basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.1_slowRate=0.5_fastRate=1_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.1 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'square'))

basename = c('estDRate_toy_agents=2000_max_time=20000_cellArea=25_fracFast=0.2_slowRate=0.5_fastRate=1_condition=constrained')
dat = rbind(dat , read_in_filelist(dir , basename , fracFree = 0.2 , cell_size = 25 , slowRate = 0.5 , fastRate = 1 , topo = 'square'))

## use the range 1000:20000 to calculate the average slope of MSD vs T, i.e. the D_eff
## note the timestep in the simulations is 0.006s, hence the constant below...
range = 1000:20000
dat %>%
  group_by(toOne, fastRate, model, fracFree , cell_size , topo) %>%
  #summarize(d_eff = mean(diff(msds[range] , 100) / 100) / 0.006) %>%
  summarize(d_eff = lm(msds[range] ~ time[range])$coefficients[2] / 0.006 / 4) %>% ## MSD ~ 4D∆t in 2 dimensions
  group_by(fastRate, model, fracFree , cell_size , topo) %>%
  mutate(d_eff = d_eff , 
         scaled_d_eff = d_eff / d_eff[1]) -> 
  d_eff_dat

unscaled_d = ggplot(d_eff_dat , aes(x = toOne , y = d_eff , col = (fracFree) , group = fracFree )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  coord_cartesian(xlim = c(0,1)) +
  ylab('d_eff (µm2/s)')
unscaled_d

scaled_d = ggplot(d_eff_dat , aes(x = toOne , y = scaled_d_eff , col = (fracFree) , group = fracFree )) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  coord_cartesian(xlim = c(0,1)) +
  ylab('d_eff (rel. tx=0)')
scaled_d

d_eff_plots = gridExtra::arrangeGrob(unscaled_d , scaled_d , nrow = 1)
#ggsave('./model/d_eff_by_fracFree.pdf' , d_eff_plots , width = length(d_eff_plots) * 4.5 , height = 3)




## plot comparision for different "fast rates"
con_d_cellsize = ggplot(d_eff_dat[d_eff_dat$model == 'constrained',] , aes(x = toOne , y = scaled_d_eff , col = (fastRate) , group = fastRate)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  #geom_smooth(span = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  #  coord_cartesian(ylim =c(2.5,6)) +
  ylab('d_eff (rel toOne = 0)') +
  ggtitle('jump-limited')
con_d_cellsize
#ggsave('./model/d_eff_by_fastRate.pdf' , con_d_cellsize , width = 4 , height = 4)


## plot comparision for different "frac free"
con_d_cellsize = ggplot(d_eff_dat[d_eff_dat$model == 'constrained',] , aes(x = toOne , y = scaled_d_eff , col = (fracFree) , group = fracFree)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  #geom_smooth(span = 0.5) +
  theme(legend.position = c(0.8,0.15)) +
  theme_classic()+
  #  coord_cartesian(ylim =c(2.5,6)) +
  ylab('d_eff (rel toOne = 0)') +
  ggtitle('jump-limited')
con_d_cellsize
#ggsave('./model/d_eff_by_fracFree.pdf' , con_d_cellsize , width = 4 , height = 4)





require(R.matlab)
require(dplyr)
require(ggplot2)

build_flux_data =  function(data , old_data = NULL , cell_size = 0){
  new_dat = NULL
  n_sets = dim(data[[1]])[2] 
  for(i in 1:n_sets){
    this_set = (data[[1]][1,i][1][[1]][[1]])
    
    to_add = data.frame(time = 1:(unlist(this_set['max.time',1,1])) , 
                        fluxxed = this_set['fluxxed',1,1],
                        model = strsplit(unlist(this_set['name',1,1]),'_')[[1]][1],
                        toTwo = as.numeric(this_set[['matrix',1,1]][1,2]) , 
                        toOne = as.numeric(this_set[['matrix',1,1]][2,1]) , 
                        cell_size = cell_size)
    new_dat = rbind(new_dat , to_add)
  }
  
  dat = rbind(old_data , new_dat)
  
  return(dat)
}


dir = '/path/to/matlab/output/'
basename = 'estFluxRate_toy_agents=2000_max_time=40000_cellSize=5_fracFast=0.05_slowRate=0.5_fastRate=1_map'

f_list = list.files(dir , full.names = T , pattern = basename)
f_list
dat = NULL
for(i in 1:length(f_list)){
  print(i)
  dat = build_flux_data(data = readMat(f_list[i]) , old_data = dat , cell_size = 5)
}

pdat = dat[dat$time %% 20 == 0 & dat$cell_size == 5,]
ff = ggplot(data = pdat[pdat$model =='free',] , aes(x = time ,y = fluxxed , col = (toTwo))) +
  geom_point(size = 0.1) +
  coord_cartesian(ylim = c(0,100000) , xlim = c(0,30000))+
  theme(legend.position = c(0.25,0.75))+
  theme_classic()+
  ggtitle('Free')
fc = ggplot(data = pdat[pdat$model == 'constrained',] , aes(x = time, y = fluxxed, col = (toTwo))) +
  geom_point(size = 0.1)+
  coord_cartesian(ylim = c(0,100000) , xlim = c(0,30000))+
  theme(legend.position = c(0.25,0.75)) +
  theme_classic()+
  ggtitle('constrained')
raw_flux_plots = gridExtra::arrangeGrob(ff , fc , nrow = 1)
plot(raw_flux_plots)
#ggsave('./model/raw_flux.pdf' , raw_flux_plots , width = 8 , height = 4)


## use the range 10000:40000 to calculate the average flux over a 100-step interval
## 0.006 comes from simulation timestep
dat= dat[dat$toTwo > 0,]
range = 10000:40000
dat %>%
  group_by(toTwo , toOne, model , cell_size) %>%
  mutate(f_eff = mean(diff(fluxxed[range] , lag = 100))/100) %>%
  summarise(f_eff = mean(f_eff) / 0.006 / 2000) %>% #2k because 2k molecules simulated
  group_by(model , cell_size)%>%
  mutate(f_eff = f_eff ,
         scaled_f_eff = f_eff / f_eff[1]) ->
  f_eff_dat
## the secretion rate in constrained , toTwo = 0 is negative over long time scales. not crazy. but it requires me to nomralize to the second value

unscaled_f = ggplot(f_eff_dat[f_eff_dat$cell_size == 5,] , aes(x = toOne , y = (f_eff) , col = model , group = model)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  #geom_smooth(alpha = 0.5 , span = 0.9)+
  theme(legend.position = c(0.8,0.2)) +
  coord_cartesian(xlim = c(0,0.25)) +
  theme_classic()+
  ylab('flux per sender-molecule per second')
unscaled_f

scaled_f = ggplot(f_eff_dat[f_eff_dat$cell_size == 5,] , aes(x = toOne , y = scaled_f_eff , col = model , group = model)) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.5) +
  #geom_smooth(alpha = 0.5 , span = 0.9)+
  theme(legend.position = c(0.8,0.2)) +
  theme_classic()+
  coord_cartesian(xlim = c(0,0.25)) +
  ylab('Flux (rel 0.005%)')
scaled_f 

f_eff_plots = gridExtra::arrangeGrob(unscaled_f , scaled_f , nrow = 1)
#ggsave('./model/f_eff_model_comparison.pdf' , f_eff_plots , width = 8 , height = 3)
 

## cell size comparison below
free_f_cellsize = ggplot(f_eff_dat[f_eff_dat$model == 'free',] , aes(x = toTwo , y = f_eff , col = factor(cell_size) , group = factor(cell_size))) +
  geom_point() +
  geom_line() +
  theme(legend.position = c(0.2,0.85)) +
  ylab('flux per sender-molecule per second') +
  ggtitle('unconstrained')
free_f_cellsize

con_f_cellsize = ggplot(f_eff_dat[f_eff_dat$model == 'constrained',] , aes(x = toTwo , y = f_eff , col = factor(cell_size) , group = factor(cell_size))) +
  geom_point() +
  geom_line() +
  theme(legend.position = c(0.2,0.85)) +
  ylab('flux per sender-molecule per second') +
  ggtitle('jump-limited')
con_f_cellsize
f_eff_plots_by_cellSize = gridExtra::arrangeGrob(free_f_cellsize , con_f_cellsize , nrow = 1)

##ggsave('./model/f_eff_model_comparison_by_cellSize.pdf' , f_eff_plots_by_cellSize , width = 8 , height = 4)

