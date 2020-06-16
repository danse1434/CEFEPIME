#-------------------------------------------------------------------------------#
# Gráficos bidimensionales -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Carga de paquetes
require(patchwork)
require(rlang)
require(tidyverse)

# Definición de directorio principal
aux_dir =  file.path(getwd(), '2_bivariado')

# Lectura de archivo funciones
source('src/6_funciones.R', encoding = 'UTF-8')

# Apertura de archivo de datos de parámetros de modelo base final
populationParameters <-
  read_csv("../4_Minimizacion/BASE_MODEL/populationParameters.txt")

popParam <- populationParameters %>% 
  select(parameter, value) %>% 
  pivot_wider(names_from = parameter, values_from = value)

# Extracción de los datos
plot1 <- extractor('Cl_pop-V1_pop', n = 900)
plot2 <- extractor('V1_pop-V2_pop', n = 900)

plot1b <- slice(plot1, which.min(`-2LL`))
# Mínimo(21.695, 18.479); t0 (20.6, 23.8)

plot2b <- slice(plot2, which.min(`-2LL`))
# Mínimo(18.479, 14.908); t0 (23.8, 13.2)


#-------------------------------------------------------------------------------#
# Gráfico 

G_plot_1 <- plot1 %>% 
  ggplot(aes(x = Cl_pop, y = V1_pop, z = LL1)) + 
    geom_contour_filled(aes(fill = stat(level)),
                        breaks = c(seq(0, 2400, by = 200))) +
  guides(fill = guide_colorsteps(barheight = unit(4, "cm"))) +
  scale_fill_viridis_d(name = 'LL') +
  xlab("Cl") + ylab(expression(V[1])) + 
  geom_point(data = popParam, aes(Cl_pop, V1_pop), shape = 4, col = 'green1', 
             inherit.aes = FALSE)  +
  geom_point(data = plot1b, shape = 8, col = 'red')


G_plot_2 <- plot2 %>% 
  ggplot(aes(x = V2_pop, y = V1_pop , z = LL1)) + 
  geom_contour_filled(aes(fill = stat(level)),
                      breaks = c(seq(0, 2400, by = 100))) +
  guides(fill = guide_colorsteps(barheight = unit(4, "cm"))) +
  scale_fill_viridis_d(name = 'LL') +
  xlab(expression(V[2])) +ylab(expression(V[1])) + 
  geom_point(data = popParam, aes(V2_pop, V1_pop), shape = 4, col = 'green1', 
             inherit.aes = FALSE)  +
  geom_point(data = plot2b, shape = 8, col = 'red')

G_plot = (G_plot_1 + G_plot_2) &
  plot_annotation(tag_levels = 'A')

# Almacenamiento
ggsave(filename = 'figures/2_perfiles_LL_BV.pdf', plot = G_plot, device = 'pdf', 
       width = 8*1.2, height = 4*1.2)




