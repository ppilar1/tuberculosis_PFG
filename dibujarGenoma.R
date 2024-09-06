library(circlize) # Para crear los gráficos circulares del genoma
library(genbankr) # Para trabajar con archivos de GenBank -> la base de datos
library(readxl) # Para leer el archivo Excel con las posiciones de las IS
library(dplyr) # Para manipular los datos

ruta_archivo <- "posiciones.xlsx" # Se guarda en la variable el Excel de las posiciones de las IS

nombres_pestanas <- excel_sheets(ruta_archivo) # Se guardan los nombres de cada una de las hojas (muestras) del Excel

lista_datos <- list() # Creamos una lista vacía para almacenar los datos de las hojas

for (pestana in nombres_pestanas) {
  lista_datos[[pestana]] <- read_excel(ruta_archivo, sheet = pestana)
} # Se lee cada hoja y se almacena su contenido en 'lista_datos'

# lista_datos$`20ID00937_S7 ` -> en caso de querer visualizar la información de una hoja en concreto


gb_file <- "~/Descargas/h37rv.gb" # Se accede al genoma completo de M. tuberculosis H37Rv 
genome_data <- readGenBank(gb_file) # Se guarda el genoma de M.tuberculosis en 'genome_data'

genome_length <- genome_data@sequence$H37Rv@length # Se guarda la longitud del genoma (~ 4,4 millones pb)

genes <- genome_data@cds # Se guardan las secuencias codificantes de genes (CDS) del genoma, que es lo que interesa


nombre.df <- names(lista_datos)[1] # Se guardan los nombres de las hojas de Excel, pero a través de la lista que se ha llenado (se podría mirar de hacer con 'nombres_pestanas)

for (nombre.df in names(lista_datos)){ # Se itera sobre cada hoja del Excel (cada elemento de la lista)
  
  gene_data <- data.frame( # Se genera un data frame
    Chromosome = rep("chr1", length(genes)), # Todas las columnas tienen el cromosoma 'chr1'
    Start = genes@ranges |> data.frame() |> select(start) |> pull(), # Columna para la posición inicial de cada gen en el genoma
    End  = genes@ranges |> data.frame() |> select(end) |> pull(),# Columna para la posición final de cada gen en el genoma
    Gene = genes$gene # Columna para los nombres de cada gen
  )
  
  gene_data <- gene_data |> filter(!is.na(Gene)) # Se eliminan las filas que no tienen valor (NA) en la columna 'Gene'
  
  gene_data<- gene_data |> slice(seq(1,nrow(gene_data),25)) # Se va a enseñar solo 1 de cada 25 genes en el gráfico

 
  gene_starts <- gene_data$Start # Se extrae la posición de inicio de cada gen que hay en el data frame
  gene_ends <- gene_data$End # Se extrae la posición final de cada gen que hay en el data frame
  gene_names <- gene_data$Gene # Se extrae el nombre de cada gen que hay en el data frame
  
  
  png(file = gsub(" ","",paste0(nombre.df,".png",sep ="")), width = 800, height = 800) # Se genera un archivo (.png) con el nombre correspondiente de cada muestra para guardar luego el gráfico circular
  
  # DIBUJO DEL GRÁFICO
  circos.par("start.degree" = 90)
  circos.initialize(factors = "chr1", xlim = c(0, genome_length),start)
  
  circos.trackPlotRegion(factors = "chr1", ylim = c(0, 1), track.height = 0.1,
                         bg.col = "lightblue", bg.border = NA) # Se genera los círculos del genoma y se le añade color
 
  data.frame.is <-lista_datos[[nombre.df]] |> 
    select(Comienzo,Final,Gen) |> 
    data.frame() # Del Excel con las posiciones se sacan las columnas 'Comienzo', 'Final' y 'Gen' para hacer otro data frame
 
  
  # Dibujar los genes como rectángulos
  # circos.rect(xleft = gene_starts,
  #             ybottom = 0, 
  #             xright = gene_ends, 
  #             ytop = 1,
  #             col = "darkgreen", 
  #             border = "white")
  
  
  gene_midpoints <- (gene_starts + gene_ends)/2 # Las etiquetas de los genes se colocarán en el medio de la extensión de estos genes
  
  # Se colocan las etiquetas de los genes en los que la IS6110 se ha insertado alrededor del círculo
  circos.text(x = gene_midpoints, y = 1.2, labels = gene_names, facing = 'clockwise', niceFacing = TRUE, cex = 1, adj = c(0, 0)) 
  circos.labels(rep("chr1",nrow(data.frame.is)), line_col = 'red', x = data.frame.is$Comienzo, labels = data.frame.is$Comienzo, cex=1)
  
  # circos.trackPlotRegion(factors = "chr1", ylim = c(0, 1), track.height = 0.8,
  #                        bg.col = "lightblue", bg.border = NA)
  
  
  #i <- 1
  # for (i in 1:nrow(data.frame.is)){
  #   circos.lines(x = c(data.frame.is[i,"Comienzo"],
  #                      data.frame.is[i,"Final"]), 
  #                y = c(1, 2.8), 
  #                col = "red", lwd = 4)
  # }
  
  # Marcas y etiquetas en el genoma
  tick_interval <- 500000  # 0.5 Mbp en pares de bases
  ticks <- seq(0, genome_length, by = tick_interval)
  ticks <- data.frame(ticks)
  
 # Mostrar las marcas y etiquetas
  circos.trackPlotRegion(factors = "chr1",
                         ylim = c(1.3, 1.4),
                         track.height = 0.05,
                         bg.col = "lightgray",
                         bg.border = NA)
  
  circos.labels(rep("chr1",nrow(ticks)), 
                line_col = 'black',
                x = ticks$ticks, 
                labels = paste0(ticks$ticks / 1e6, " Mbp"),
                cex=2)
  
  text(0, 0, nombre.df, cex = 2) # Añade el nombre de la muestra al centro del gráfico (dentro del círculo)
  
  circos.clear() # Reestablecimiento del gráfico y de la función
  dev.off() # Se cierra el gráfico, guardándolo en el archivo del nombre de la muestra al que corresponde
}


