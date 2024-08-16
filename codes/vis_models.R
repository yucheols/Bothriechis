##### visualize model outputs

# clear working environments
rm(list = ls(all.names = T))
gc()

# load packages
library(terra)
library(rasterVis)
library(ggplot2)
library(viridis)

##### continuous model ----------
cont <- rast(list.files(path = 'output_models_preds', pattern = '.tif$', full.names = T))
cont <- terra::subset(cont, c(2,4))

names(cont) = c('Green', 'Yellow')
print(cont)

## plot
gplot(cont) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  scale_fill_gradientn(colors = inferno(1000),
                       na.value = NA,
                       name = 'Suitability',
                       breaks = c(0.1, 0.9),
                       labels = c('Low: 0.1', 'High: 0.9')) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.title = element_text(size = 14, margin = margin(b = 12), face = 'bold'),
        legend.text = element_text(size = 12))

# save
ggsave('plots/Bothriechis_cont.png', width = 25, height = 10, dpi = 800, units = 'cm')


##### binary model ----------
bin <- rast(list.files(path = 'output_models_preds', pattern = '.tif$', full.names = T))
bin <- terra::subset(bin, c(1,3))

names(bin) = c('Green', 'Yellow')
print(bin)

## plot
gplot(bin) +
  geom_tile(aes(fill = value)) +
  coord_equal() +
  facet_wrap(~ variable) +
  xlab('Longitude (°)') + ylab('Latitude (°)') +
  scale_fill_gradientn(colors = rev(terrain.colors(1000)),
                       na.value = NA) +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 12)),
        legend.position = 'none')

# save
ggsave('plots/Bothriechis_bin.png', width = 25, height = 10, dpi = 800, units = 'cm')


