setwd("PATH DONDE TENGAS ESTE CODIGO Y LAS SECUENCIAS") # Se establece el directorio de trabajo en la ruta especificada en la que se quiera trabajar

ficheros.fastq<-list.files("../../Raw/",pattern="fastq",full.names=T) # Se guardan los ficheros FASTQ (".fastq") en una lista
is6110.inicio<-toupper("gagagtctccggactcaccggggcggttca") # Se guarda la IS6110 correspondiente en mayúsculas
longitud.target.inicio<-nchar(is6110.inicio) # Se guarda la longitud de la IS6110 a través del contaje del número de nucléotidos que tiene

## Docker ShortRead
library(ShortRead)

# Se quiere probar solamente con unos de los ficheros para hacerlo más simple
ficheros.fastq<-c("HMS14032_lib66_R1.fastq","HMS2742_lib102_R1.fastq") # Se selecciona dos archivos en el vector (estos archivos pueden variar)

lista.data.frames<-lapply(ficheros.fastq,FUN=function(fichero){
  fq <- readFastq(fichero) # Se lee el archivo FASTQ
  fq.text<-sread(fq) # Se extraen las secuencias de ADN
  buscar.patron.inicio<-start(vmatchPattern(pattern = is6110.inicio, subject = fq.text, fixed=TRUE)) # Se encuentra el inicio de las coincidencias de la cadena IS en cada uno de los reads
  detectar.reads.inicio<-sapply(buscar.patron.inicio,FUN=function(x){ length(x)!=0}) # Se detectan las secuencias que contienen el patrón de búsqueda
  names(buscar.patron.inicio)<-1:length(fq.text) # Se asignan nombres a las posiciones de 'buscar.patron.inicio'
  indice.fastq.detectadas.inicio<-as.numeric(paste0(names(buscar.patron.inicio[detectar.reads.inicio])))  # Se obtienen los índices de las secuencias donde se ha detectado el patrón de búsqueda
  fastq.detectadas.inicio<-fq.text[indice.fastq.detectadas.inicio]  # Se extran las secuencias en donde se ha detectado el patrón de búsqueda a través de los índices
  posicion.secuencia.inicio.en.fastq<-buscar.patron.inicio[indice.fastq.detectadas.inicio] # Se obtienen las posiciones de donde se ha encontrado la coincidencia del patrón de búsqueda en las secuencias
  list(fastq.detectectadas.inicio=fastq.detectadas.inicio,posicion.secuencia.inicio.en.fastq=posicion.secuencia.inicio.en.fastq) # Se devuelve una lista con las secuencias en donde se ha detectado el patrón de búsuqeda y las posiciones correspondientes donde se encuentra el patrón
})

names(lista.data.frames)<-basename(ficheros.fastq) # Cada elemento de 'lista.data.frames' tiene el nombre de su fichero FASTQ correspondiente
save(lista.data.frames,file="lista.data.frames.final.Rdata") # Se accede al elemento de 'lista.data.frames' con ese nombre concreto


lista.resultados<-lapply(lista.data.frames,FUN=function(df){
  num.secuencias.df<-length(df[[1]]) # Se cuenta el número de lecturas que tiene 'df'
  
  # Se obtiene qué hay desde el principio de la IS6110 hasta el final de la secuencia
  target.inicio.final<-lapply(1:num.secuencias.df,FUN=function(indice.fastq){
    secuencia<-df[[1]][[indice.fastq]]
    indicador.comienzo.is<-df[[2]][[indice.fastq]]
    longitud.secuencia<-length(secuencia)
    toString(secuencia[indicador.comienzo.is:longitud.secuencia])
  })
  
  # Se extraen los 6 nucleotidos que hay justo después  de la IS6110 
  post.6.nucleotidos<-lapply(1:num.secuencias.df,FUN=function(indice.fastq){
    secuencia<-df[[1]][[indice.fastq]]
    indicador.comienzo.is<-df[[2]][[indice.fastq]]
    longitud.secuencia<-length(secuencia)
    ifelse(longitud.secuencia-(indicador.comienzo.is+longitud.target.inicio)>=5,
           toString(secuencia[(indicador.comienzo.is+longitud.target.inicio):(indicador.comienzo.is+longitud.target.inicio+5)]),
           NA) # Si la IS6110 coincide con el final de la secuencia y no hay 6 nucléotidos posteriores, se pone NA
  })
  
  # Se extrae la parte de antes de la IS6110
  pre.target<-lapply(1:num.secuencias.df,FUN=function(indice.fastq){
    secuencia<-df[[1]][[indice.fastq]]
    indicador.comienzo.is<-df[[2]][[indice.fastq]]
    longitud.secuencia<-length(secuencia)
    ifelse(indicador.comienzo.is==1,NA,
           toString(secuencia[1:(indicador.comienzo.is-1)]))
  })
  
  list(target.inicio.final=target.inicio.final,
       post.6.nucleotidos=post.6.nucleotidos,
       pre.target=pre.target)
  # Se genera una lista con los elementos especificados -> hay 3 elementos dentro de la lista cada uno correspondiente con las 3 zonas identificadas
  
})

names(lista.resultados)<-basename(ficheros.fastq) # Cada elemento de 'lista.resultados' tiene el nombre de su fichero FASTQ correspondiente
save(lista.resultados,file="lista.resultados.final.is.rdata") # Se guarda la variable 'lista.resultados' en un archivo con extensión '.Rdata'

## Docker tidyverse
library(tidyverse)
load("lista.resultados.final.is.rdata") # Se carga el archivo '.Rdata' que se acaba de guardar

# Se cuenta cuántas veces aparece cada cadena de 6 nucleótidos posterior a la IS6110 en los reads
resultado<-lista.resultados[[1]] # Se selecciona solamente el primer elemento de 'lista.resultados'
lista.post.6.nucleotidos<-lapply(lista.resultados,FUN=function(resultado){
  data.frame.6.nucleotidos<-do.call("rbind.data.frame",resultado$post.6.nucleotidos) # Se juntan los datos en un solo data frame
  colnames(data.frame.6.nucleotidos)<-"fragmento" # Se renombra la columna del data frame
  data.frame.6.nucleotidos
  sort(table(data.frame.6.nucleotidos$fragmento),decreasing = T) %>% # Se cuenta la frecuencia de aparición de los sextetos identificados y se ordena de forma descendente
    data.frame() # El resultado es de tipo data frame 
})

names(lista.post.6.nucleotidos)<-basename(ficheros.fastq) # Cada elemento de 'lista.post.6.nucleotidos' tiene el nombre de su fichero FASTQ correspondiente

# Se ordena la parte que hay posterior a la IS6110 
lista.target.inicio.final<-lapply(lista.resultados,FUN=function(resultado){
  data.frame.target.inicio.final<-do.call("rbind.data.frame",resultado$target.inicio.final) # Se juntan los datos en un solo data frame
  colnames(data.frame.target.inicio.final)<-"secuencia.target.inicio.final" # Se renombra la columna del data frame
  longitud.cadena<-stringr::str_length(data.frame.target.inicio.final$secuencia.target.inicio.final) # Se calcula la longitud de los fragmentos identificados
  data.frame.target.inicio.final<-data.frame.target.inicio.final[order(longitud.cadena,decreasing = T),] # Se ordenan los resultados en función de la longitud anterior y de forma descendente
  data.frame.target.inicio.final # El resultado es de tipo data frame 
})

names(lista.target.inicio.final)<-basename(ficheros.fastq) # Cada elemento de 'lista.target.inicio.final' tiene el nombre de su fichero FASTQ correspondiente

# Se ordena la parte que hay desde el comienzo IS6110 hasta el final de toda la secuencia
lista.pre.target<-lapply(lista.resultados,FUN=function(resultado){
  data.frame.pre.target<-do.call("rbind.data.frame",resultado$pre.target) # Se juntan los datos de 'pre.target' en un solo data frame
  data.frame.target.inicio.final<-do.call("rbind.data.frame",resultado$target.inicio.final) # Se juntan los datos de 'target.inicio.final' en un solo data frame
  colnames(data.frame.pre.target)<-"secuencia.pre.target" # Se renombra la columna 'pre.target' del data frame
  longitud.cadena<-stringr::str_length(data.frame.pre.target$data.frame.target.inicio.final) # Se calcula la longitud de cada fragmento de la columna renombrada
  data.frame.pre.target<-data.frame.pre.target[order(longitud.cadena,decreasing = T),] # Se ordena el marco de datos 'data.frame.pre.target' según la longitud de cadena (que es la longitud de la parte correspondiente al pre.target) en orden descendente. Esto significa que las filas con cadenas más largas aparecerán primero en el marco de datos resultante
})

names(lista.pre.target)<-basename(ficheros.fastq) # Cada elemento de 'lista.pre.target' tiene el nombre de su fichero FASTQ correspondiente

lista.resultados.en.data.frame<-lapply(lista.resultados,FUN=function(muestra){
  df.post<-do.call("rbind.data.frame",muestra[[1]]) # Se juntan los datos de post.target en un solo data frame que corresponden con la primera columna 
  df.seis<-do.call("rbind.data.frame",muestra[[2]]) # Se juntan los datos de los 6NT posteriores a la IS en un solo data frame que corresponden con la segunda columna
  df.pre<-do.call("rbind.data.frame",muestra[[3]]) # Se juntan los datos de pre.target en un solo data frame que corresponden con la tercera columna
  length.pre<-stringr::str_length(df.post[,1]) # Se calcula la longitud de los fragmentos de tipo post.target
  df.resultado.muestra<-data.frame(df.post,length.pre,df.seis,df.pre) # Se genera un nuevo data frame que combina los tres anteriores creados junto con las longitudes
  colnames(df.resultado.muestra)<-c("Target_Final","Length","Seis.Nucleotidos","Pre") # Se les da nombre a las columnas del nuevo data frame
  df.resultado.muestra<-df.resultado.muestra[order(df.resultado.muestra$Length,decreasing = T),] # Se ordena en función de los valores de la columna 'Length' de forma descendente
  df.resultado.muestra # El resultado es de tipo data frame 
})

names(lista.resultados.en.data.frame)<-basename(ficheros.fastq) # Cada elemento de 'lista.resultados.en.data.frame' tiene el nombre de su fichero FASTQ correspondiente

# En 'lista.resultados.en.data.frame' tenemos un data frame con cuatro columnas: la parte post, la longitud de la parte previa a la IS6110, los seis nucleótidos posteriores a la IS y la parte pre

# Se limpian los datos
lista.resultados.filtrados<-lapply(lista.resultados.en.data.frame,FUN=function(muestra){
  tabla.seis.nucleotidos<-
    muestra %>% 
    dplyr::count(Seis.Nucleotidos) %>% # Se cuenta el número de veces que cada sexteto posterior a la IS aparece en 'muestra'
    dplyr::filter(!is.na(Seis.Nucleotidos)&(n>1)) %>%  # Se escogen aquellos resultados que no son 'NA' y que aparecen más de una vez
    arrange(desc(n)) # Se ordena de forma descendente
  
  lapply(tabla.seis.nucleotidos$Seis.Nucleotidos,FUN=function(cadena){
    muestra %>% 
      dplyr::filter(Seis.Nucleotidos%in%cadena) %>% # Se filtran las filas de 'muestra' que contienen el sexteto seleccionado en 'cadena'
      arrange(desc(Length)) %>%  # Se ordenan esas filas en función de la longitud que tienen y de forma descendente
      dplyr::slice(1:3) # Se eligen las tres primeras filas de la lista ordenada
  })
})

names(lista.resultados.filtrados)<-basename(ficheros.fastq) # Cada elemento de 'lista.resultados.filtrados' tiene el nombre de su fichero FASTQ correspondiente
save(lista.post.6.nucleotidos,lista.target.inicio.final,lista.pre.target,file="resultados.finales.final.is.Rdata") # Se guardan 'lista.post.6.nucleotidos', 'lista.target.inicio.final' y 'lista.pre.target.final' en un archivo con extensión '.Rdata'


load("resultados.finales.final.is.Rdata")
ficheros.fastq<-list.files("../../Raw/",pattern="fastq",full.names=T)

# Creación de la estructura de los futuros archivos Excel
nficheros<-length(ficheros.fastq) # Se cuenta el número de ficheros FASTQ
nombres<-basename(ficheros.fastq) # Se guardan los nombres de los ficheros sin la extensión
nexcel<-trunc(161/10)

lista.para.excel<-lapply(1:nexcel,FUN=function(x){
  v<-1+((x*10)-10)
  seq(v,length.out = 10,by = 1)
}) # Los archivos Excel tendrán datos de hasta 10 ficheros FASTQ

lista.para.excel[[17]]<-161
names(lista.para.excel)<-1:17

lista.para.excel
indice.excel<-13
indice.fichero<-lista.para.excel[[indice.excel]]
documento<-indice.fichero[[1]]

library(xlsx)

names(lista.pre.target.inicio)[121:130]

lapply(1:length(lista.para.excel),FUN=function(indice.excel){
  lapply(lista.para.excel[[indice.excel]],FUN=function(indice.fichero){
    lapply(indice.fichero,FUN=function(documento){
       xlsx::write.xlsx(x = lista.post.6.nucleotidos[[documento]],file = paste0("post6.nucleotidos",indice.excel,".xlsx"),append=T,sheetName = names(lista.post.6.nucleotidos)[documento]) # Se crea un Excel en donde sale la información de 'lista.post.6.nucleotidos' para cada FASTQ
       xlsx::write.xlsx(x = lista.target.inicio.final[[documento]],file = paste0("excel.target.inicio.final",indice.excel,".xlsx"),append=T,sheetName = names(lista.target.inicio.final)[documento]) # Se crea un Excel en donde sale la información de 'excel.target.inicio.final' para cada FASTQ
       xlsx::write.xlsx(x = lista.pre.target[[documento]],file = paste0("excel.pre.target",indice.excel,".xlsx"),append=T,sheetName = names(lista.pre.target)[documento]) # Se crea un Excel en donde sale la información de 'excel.pre.target' para cada FASTQ 
    })  
    })
  })

# Escribimos la parte `target.inicio.final` en ficheros FASTA
ficheros.fastq<-list.files("../../Raw/",pattern="fastq",full.names=T)
lapply(1:length(lista.target.inicio.final),FUN=function(cepa){
  nombre.fichero<-paste0(names(lista.target.inicio.final)[cepa],".fasta")
  lapply(1:length(lista.target.inicio.final[[cepa]]),FUN=function(secuencia){
    seqinr::write.fasta(sequences = as.character(paste0(lista.target.inicio.final[[cepa]]))[secuencia],
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





