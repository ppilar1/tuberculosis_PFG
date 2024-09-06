setwd("PATH DONDE TENGAS ESTE CODIGO Y LAS SECUENCIAS") # Se establece el directorio de trabajo en la ruta especificada en la que se quiera trabajar

# Se quiere detectar qué ficheros FASTQ hay
ficheros.fastq<-list.files(path = ".",pattern=".fastq$",full.names=T) # Se guardan los ficheros FASTQ (".fastq") en una lista
is6110.inicio<-toupper("tgaaccgccccggcatgtccggagactccagtt") # Se guarda la IS6110 correspondiente en mayúsculas
longitud.target.inicio<-nchar(is6110.inicio) # Se guarda la longitud de la IS6110 a través del contaje del número de nucléotidos que tiene

## Docker ShortRead
library(ShortRead)

# Se quiere probar solamente con unos de los ficheros para hacerlo más simple
ficheros.fastq <- ficheros.fastq[1:2] # Se selecciona los primeros dos archivos en el vector
fichero <- ficheros.fastq[1] # Se asigna el primer archivo de tipo FASTQ a la variable 'fichero'

lista.data.frames<-lapply(ficheros.fastq,FUN=function(fichero){
  fq <- readFastq(fichero) # Se lee el archivo FASTQ
  fq.text<-sread(fq) # Se extraen las secuencias de ADN
  class(fq.text) # Se verifica que la 'fq.text' es de la clase 'ShortRead'
  buscar.patron.inicio<-start(vmatchPattern(pattern = is6110.inicio, subject = fq.text, fixed=TRUE)) # Se encuentra el inicio de las coincidencias de la cadena IS en cada uno de los reads
  # class(buscar.patron.inicio)
  detectar.reads.inicio<-sapply(buscar.patron.inicio,FUN=function(x){ length(x)!=0}) # Se detectan las secuencias que contienen el patrón de búsqueda
  names(buscar.patron.inicio)<-1:length(fq.text) # Se asignan nombres a las posiciones de 'buscar.patron.inicio'
  indice.fastq.detectadas.inicio<-as.numeric(paste0(names(buscar.patron.inicio[detectar.reads.inicio]))) # Se obtienen los índices de las secuencias donde se ha detectado el patrón de búsqueda
  fastq.detectadas.inicio<-fq.text[indice.fastq.detectadas.inicio] # Se extran las secuencias en donde se ha detectado el patrón de búsqueda a través de los índices
  posicion.secuencia.inicio.en.fastq<-buscar.patron.inicio[indice.fastq.detectadas.inicio] # Se obtienen las posiciones de donde se ha encontrado la coincidencia del patrón de búsqueda en las secuencias
  list(fastq.detectectadas.inicio=fastq.detectadas.inicio,posicion.secuencia.inicio.en.fastq=posicion.secuencia.inicio.en.fastq) # Se devuelve una lista con las secuencias en donde se ha detectado el patrón de búsuqeda y las posiciones correspondientes donde se encuentra el patrón
})

names(lista.data.frames)<-basename(ficheros.fastq) # Cada elemento de 'lista.data.frames' tiene el nombre de su fichero FASTQ correspondiente
lista.data.frames$`20ID00931_S1_R1_001.fastq` # Esta línea y la de abajo hacen lo mismo: se accede al elemento de 'lista.data.frames' con ese nombre concreto
#lista.data.frames[['20ID00931_S1_R1_001.fastq']]

#save.image(file = "environment_150324.Rdata") Se guardará todo el entorno en el que se esté trabajando en ese momento en un archivo con extensión '.Rdata'
#load("environment_150324.Rdata") Se carga el entorno de trabajo que se acaba de guardar
save(lista.data.frames,file="lista.data.frames.Rdata") # Se guarda la variable 'lista.data.frames' en un archivo con extensión '.Rdata'

# Se quiere probar solamente con uno de los objetos de 'lista.data.frames' para hacerlo más simple
df <- lista.data.frames[[1]] # Se quiere trabajar solamente con el primer elemento
length(lista.data.frames) # Se calcula el número de elementos que contiene la lista
indice.fastq = 1
df[[2]][37] # Se accede al valor número 37  del segundo componente del primer elemento de lista

lista.resultados<-lapply(lista.data.frames,FUN=function(df){
  num.secuencias.df<-length(df[["fastq.detectectadas.inicio"]]) # Se cuenta el número de lecturas que tiene 'df'
  
  # Se obtiene qué hay desde el principio de la secuencia hasta el final de la IS6110
  pre.target.inicio<-lapply(1:num.secuencias.df,FUN=function(indice.fastq){
    secuencia<-df[["fastq.detectectadas.inicio"]][[indice.fastq]]
    indicador.comienzo.is<-df[["posicion.secuencia.inicio.en.fastq"]][[indice.fastq]]
    longitud.secuencia<-length(secuencia)
    toString(secuencia[1:(indicador.comienzo.is+longitud.target.inicio-1)])
  })
    
  # Se extraen los 6 nucleotidos que hay justo delante de la IS6110
  pre.6.nucleotidos<-lapply(1:num.secuencias.df,FUN=function(indice.fastq){
    secuencia<-df[["fastq.detectectadas.inicio"]][[indice.fastq]]
    indicador.comienzo.is<-df[["posicion.secuencia.inicio.en.fastq"]][[indice.fastq]]
    longitud.secuencia<-length(secuencia)
    ifelse(indicador.comienzo.is>6,
           toString(secuencia[(indicador.comienzo.is-6):(indicador.comienzo.is-1)]),
           NA) # Si la IS6110 coincide con el principio de la secuencia y no hay 6 nucléotidos anteriores, se pone NA
  })
  
  # Se extrae la parte de después de la IS6110
  post.target.inicio<-lapply(1:num.secuencias.df,FUN=function(indice.fastq){
    secuencia<-df[["fastq.detectectadas.inicio"]][[indice.fastq]]
    indicador.comienzo.is<-df[["posicion.secuencia.inicio.en.fastq"]][[indice.fastq]]
    longitud.secuencia<-length(secuencia)
    ifelse(indicador.comienzo.is+longitud.target.inicio>longitud.secuencia,
           NA,
           toString(secuencia[(indicador.comienzo.is+longitud.target.inicio):(longitud.secuencia)]))
  })
  
  list(pre.target.inicio=pre.target.inicio,
       pre.6.nucleotidos=pre.6.nucleotidos,
       post.target.inicio=post.target.inicio)
  # Se genera una lista con los elementos especificados -> hay 3 elementos dentro de la lista cada uno correspondiente con las 3 zonas identificadas
  
})

names(lista.resultados)<-basename(ficheros.fastq) # Cada elemento de 'lista.resultados' tiene el nombre de su fichero FASTQ correspondiente
save(lista.resultados,file="lista.resultados.rdata") # Se guarda la variable 'lista.resultados' en un archivo con extensión '.Rdata'

## Docker tidyverse
library(tidyverse)
load("lista.resultados.rdata") # Se carga el archivo '.Rdata' que se acaba de guardar

# Se cuenta cuántas veces aparece cada cadena de 6 nucleótidos previa a la IS6110 en los reads
resultado <- lista.resultados[[1]] # Se selecciona solamente el primer elemento de 'lista.resultados'
lista.pre.6.nucleotidos<-lapply(lista.resultados,FUN=function(resultado){
  data.frame.6.nucleotidos<-do.call("rbind.data.frame",resultado$pre.6.nucleotidos) # Se juntan los datos en un solo data frame
  colnames(data.frame.6.nucleotidos)<-"fragmento" # Se renombra la columna del data frame
  data.frame.6.nucleotidos
  sort(table(data.frame.6.nucleotidos$fragmento),decreasing = T) %>% # Se cuenta la frecuencia de aparición de los sextetos identificados y se ordena de forma descendente
    data.frame() # El resultado es de tipo data frame 
})

names(lista.pre.6.nucleotidos)<-basename(ficheros.fastq) # Cada elemento de 'lista.pre.6.nucleotidos' tiene el nombre de su fichero FASTQ correspondiente
#resultado 

# Se ordena la parte que hay previa a la IS6110 hasta el final de la IS según su longitud
lista.pre.target.inicio<-lapply(lista.resultados,FUN=function(resultado){
  data.frame.pre.target<-do.call("rbind.data.frame",resultado$pre.target.inicio) # Se juntan los datos en un solo data frame
  colnames(data.frame.pre.target)<-"secuencia.pre.target" # Se renombra la columna del data frame
  longitud.cadena<-stringr::str_length(data.frame.pre.target$secuencia.pre.target) # Se calcula la longitud de los fragmentos identificados
  data.frame.pre.target<-data.frame.pre.target[order(longitud.cadena,decreasing = T),] # Se ordenan los resultados en función de la longitud anterior y de forma descendente
  data.frame.pre.target # El resultado es de tipo data frame 
})

names(lista.pre.target.inicio)<-basename(ficheros.fastq) # Cada elemento de 'lista.pre.target.inicio' tiene el nombre de su fichero FASTQ correspondiente

# Se ordena la parte que hay posterior a la IS6110 
lista.post.target<-lapply(lista.resultados,FUN=function(resultado){
  data.frame.pre.target<-do.call("rbind.data.frame",resultado$pre.target.inicio) # Se juntan los datos de 'pre.target.inicio' en un solo data frame
  data.frame.post.target<-do.call("rbind.data.frame",resultado$post.target.inicio) # Se juntan los datos de 'post.target.inicio' en un solo data frame
  colnames(data.frame.pre.target)<-"secuencia.pre.target" # Se renombra la columna 'pre.target.inicio' del data frame
  longitud.cadena<-stringr::str_length(data.frame.pre.target$secuencia.pre.target) # Se calcula la longitud de cada fragmento de la columna renombrada
  data.frame.post.target<-data.frame.post.target[order(longitud.cadena,decreasing = T),] # Se ordena el marco de datos 'data.frame.post.target' según la longitud de cadena (que es la longitud de la parte correspondiente al pre.target) en orden descendente. Esto significa que las filas con cadenas más largas aparecerán primero en el marco de datos resultante
})

names(lista.post.target)<-basename(ficheros.fastq) # Cada elemento de 'lista.post.target' tiene el nombre de su fichero FASTQ correspondiente

lista.resultados.en.data.frame<-lapply(lista.resultados,FUN=function(muestra){
  df.pre<-do.call("rbind.data.frame",muestra[[1]]) # Se juntan los datos de pre.target en un solo data frame que corresponden con la primera columna
  df.seis<-do.call("rbind.data.frame",muestra[[2]]) # Se juntan los datos de los 6NT previos a la IS en un solo data frame que corresponden con la segunda columna
  df.post<-do.call("rbind.data.frame",muestra[[3]]) # Se juntan los datos de post.target en un solo data frame que corresponden con la tercera columna
  length.pre<-stringr::str_length(df.pre[,1]) # Se calcula la longitud de los fragmentos de tipo pre.target
  df.resultado.muestra<-data.frame(df.pre,length.pre,df.seis,df.post) # Se genera un nuevo data frame que combina los tres anteriores creados junto con las longitudes
  colnames(df.resultado.muestra)<-c("Pre","Length","Seis.Nucleotidos","Post") # Se les da nombre a las columnas del nuevo data frame
  df.resultado.muestra<-df.resultado.muestra[order(df.resultado.muestra$Length,decreasing = T),] # Se ordena en función de los valores de la columna 'Length' de forma descendente
  #list(df.resultado.muestra)
  df.resultado.muestra # El resultado es de tipo data frame 
})

names(lista.resultados.en.data.frame)<-basename(ficheros.fastq) # Cada elemento de 'lista.resultados.en.data.frame' tiene el nombre de su fichero FASTQ correspondiente

# En 'lista.resultados.en.data.frame' tenemos un data frame con cuatro columnas: la parte pre, la longitud de la parte previa a la IS6110, los seis nucleótidos previos a la IS y la parte post

# Se limpian los datos
muestra <- lista.resultados.en.data.frame[[1]] # Se selecciona solamente el primer data frame de 'lista.resultados.en.data.frame'

lista.resultados.filtrados<-lapply(lista.resultados.en.data.frame,FUN=function(muestra){
  tabla.seis.nucleotidos<-
    muestra %>% 
    dplyr::count(Seis.Nucleotidos) %>% # Se cuenta el número de veces que cada sexteto previo a la IS aparece en 'muestra'
    dplyr::filter(!is.na(Seis.Nucleotidos)&(n>1)) %>% # Se escogen aquellos resultados que no son 'NA' y que aparecen más de una vez
    arrange(desc(n)) # Se ordena de forma descendente
  
  cadena <- tabla.seis.nucleotidos$Seis.Nucleotidos[1] # Se selecciona el primer sexteto más frecuente
  lapply(tabla.seis.nucleotidos$Seis.Nucleotidos,FUN=function(cadena){
    muestra %>%
      dplyr::filter(Seis.Nucleotidos%in%cadena) %>% # Se filtran las filas de 'muestra' que contienen el sexteto seleccionado en 'cadena'
      arrange(desc(Length)) %>% # Se ordenan esas filas en función de la longitud que tienen y de forma descendente
      dplyr::slice(1:3) # Se eligen las tres primeras filas de la lista ordenada
  })
})

names(lista.resultados.filtrados)<-basename(ficheros.fastq) # Cada elemento de 'lista.resultados.filtrados' tiene el nombre de su fichero FASTQ correspondiente

lista.resultados.filtrados[[1]][[1]] |> # Se selecciona el primer resultado que se ha filtrado en el primer fichero
save(lista.pre.6.nucleotidos,lista.pre.target.inicio,lista.post.target,file="resultados.finales.Rdata") # Se guardan 'lista.pre.6.nucleotidos', 'lista.pre.target.inicio' y 'lista.post.target' en un archivo con extensión '.Rdata'


# Creación de la estructura de los futuros archivos Excel
nficheros<-length(ficheros.fastq) # Se cuenta el número de ficheros FASTQ
nombres<-basename(ficheros.fastq) # Se guardan los nombres de los ficheros sin la extensión
nexcel<-trunc(nficheros/10)
nexcel<-1

lista.para.excel<-lapply(1:nexcel,FUN=function(x){
  v<-1+((x*10)-10)
  seq(v,length.out = 10,by = 1)
}) # Los archivos Excel tendrán datos de hasta 10 ficheros FASTQ

lista.para.excel[[17]]<-161
names(lista.para.excel)<-1

lista.para.excel
indice.excel<-13
indice.fichero<-lista.para.excel[[indice.excel]]
documento<-indice.fichero[[1]]

library(xlsx)

names(lista.pre.target.inicio)[121:130]

lapply(1:length(lista.para.excel),FUN=function(indice.excel){
  lapply(lista.para.excel[[indice.excel]],FUN=function(indice.fichero){
    lapply(indice.fichero,FUN=function(documento){
      xlsx::write.xlsx(x = lista.pre.6.nucleotidos[[documento]],file = paste0("pre6.nucleotidos",indice.excel,".xlsx"),append=T,sheetName = names(lista.pre.6.nucleotidos)[documento]) # Se crea un Excel en donde sale la información de 'lista.pre.6.nucleotidos' para cada FASTQ
       xlsx::write.xlsx(x = lista.pre.target.inicio[[documento]],file = paste0("excel.pre.target.inicio.",indice.excel,".xlsx"),append=T,sheetName = names(lista.pre.target.inicio)[documento]) # Se crea un Excel en donde sale la información de 'lista.pre.target.inicio' para cada FASTQ
       xlsx::write.xlsx(x = lista.post.target[[documento]],file = paste0("excel.post.target",indice.excel,".xlsx"),append=T,sheetName = names(lista.post.target)[documento]) # Se crea un Excel en donde sale la información de 'lista.post.target' para cada FASTQ 
    })  
  })
})

load("resultados.finales.Rdata") # Se cargan los datos de 'resultados.finales.Rdata' en el entorno de trabajo

# Escribimos la parte `pre.target` en ficheros FASTA
ficheros.fastq<-list.files("../../Raw/",pattern="fastq",full.names=T)
lapply(1:length(lista.pre.target.inicio),FUN=function(cepa){
  nombre.fichero<-paste0(names(lista.pre.target.inicio)[cepa],".fasta")
  lapply(1:length(lista.pre.target.inicio[[cepa]]),FUN=function(secuencia){
    seqinr::write.fasta(sequences = as.character(paste0(lista.pre.target.inicio[[cepa]]))[secuencia],
                        names=secuencia,
                        file.out = nombre.fichero,open = "a")
  })
})

ficheros.fasta<-list.files(pattern = "*.fasta",full.names = T,recursive = T,include.dirs = T) # Se escogen todos los archivos '.fasta'
nombresalida<-paste0(list.files(pattern="fasta"),".blastn") # Se crea una lista con los nombres de salida de los archivos de BLASTn. Ahora los archivos que pasen por el BLASTn tendrán la extensión '.blastn'
paste0("blastn -query ",ficheros.fasta," -db h37rv -out ", paste(nombresalida), " -outfmt 6 ", collapse = ";") # Se aplica el BLASTn a los '.fasta' utilizando la base de datos 'h37rv'

library(rtracklayer)

# Anotación de los archivos resultantes de BLASTn
nombres.columnas<-read.delim("encabezado2.txt",header=T,sep=",") # Se lee el archivo que contiene el nombre de las columnas de la anotación
ficheros.blastn<-list.files(pattern="blastn",full.names=T,recursive = T,include.dirs=T) # Se guardan las rutas completas de los archivos BLASTn
ficheros.fasta<-list.files(pattern = "*.fasta$",full.names = T,recursive = T,include.dirs = T,ignore.case = F) # Se guardan las rutas completas de los archivos FASTA
anotacion.tb<-readGFF("Mycobacterium_tuberculosis_h37rv.ASM19595v2.46.chromosome.Chromosome.gff3") # Se lee el archivo GFF correspondiente
anotacion.tb.mRNA<-anotacion.tb %>% data.frame %>% dplyr::filter(type=="mRNA") %>% dplyr::select(Parent,transcript_id,Name) # Se filtran las anotaciones para incluir solo las de tipo 'mRNA'
anotacion.tb.genes<-anotacion.tb %>% data.frame %>% dplyr::filter(type=="gene") %>% dplyr::select(ID,start,end,strand,Name,gene_id,description) # Se filtran las anotaciones para incluir solo las de tipo 'gene'
anotacion.df<-cbind(anotacion.tb.genes,subject.acc.ver=anotacion.tb.mRNA[,2]) # Combinación de data frames
anotacion.df$subject.acc.ver<-as.character(paste0(anotacion.df$subject.acc.ver))
nombres<-gsub(pattern = "fastq.gz.fasta",x = basename(ficheros.fasta),replacement = "") # Se elimina la parte '.fastq.gz.fasta' de los archivos FASTA

lapply(1:length(ficheros.fasta),FUN=function(indice){
file.blast<-read.delim(ficheros.blastn[indice],header = F)
file.fasta<-read.delim(ficheros.fasta[indice],header=F)
colnames(file.blast)<-colnames(nombres.columnas)
file.fasta<-file.fasta[seq(from=2,to=nrow(leer.fasta),by=2),]
file.fasta<-as.character(paste0(file.fasta))
df.fasta<-data.frame(cadena=file.fasta,query.acc.ver=1:length(file.fasta))
df.resultado<-left_join(df.fasta, file.blast,by="query.acc.ver")
df.resultado$subject.acc.ver<-as.character(paste0(df.resultado$subject.acc.ver))
df.resultado<-left_join(df.resultado,anotacion.df,by="subject.acc.ver")
df.resultado<-df.resultado %>% dplyr::select(2,1,3:20)
writexl::write_xlsx(df.resultado,paste0(nombres[indice],"xlsx"))
}) # Se procesa cada archivo FASTA y BLASTn, se unen a través de la anotación del genoma y se guarda el resultado en un archivo de tipo Excel
