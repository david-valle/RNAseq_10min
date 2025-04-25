# RNA-seq en 10 minutos (más o menos)

## Contacto

***David Valle-García***

[*Laboratorio de Epigenómica del Envejecimiento*](https://epiaginglab.org/)

[*Centro de Investigación sobre Envejecimiento*](https://cinvestav.mx/sur/conocenos/departamentos/cie)

*Cinvestav Sede Sur*

david.valle -at- cinvestav.mx

## Requerimientos

En este documento están los pasos básicos para realizar un análisis rápido de RNA-seq. Para el ejericicio, es necesario descargar las carpetas de [este link](https://drive.google.com/drive/folders/1Dg9OvEoRqitZiIHejZQLZgopR7ZhysPl). Guárdalas todas en el directorio donde vayas a realizar el ejercicio. Todos los comandos asumen que tienes esas carpetas con esos nombres en el directorio donde los estás corriendo.

Adicionalmente, necesitas instalar los siguientes paquetes.

En la terminal (vía [conda](https://conda-forge.org/)):

-   salmon

En R (studio):

-   edgeR

-   ComplexHeatmap

## Introducción

El análisis de RNA-seq es una metodología cada vez más accesible que nos permite descubrir una gran cantidad de patrones biológicos. Dada su naturaleza, es difícil hacer un método global que abarque todos los casos, pero en este ejercicio analizaremos el tipo de análisis más simple: comparar entre dos condiciones biológicas.

## Los datos

Los datos que analizaremos provienen del dataset [GSE168137](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168137). Son datos pareados de RNA-seq de ratones silvestres (tanto hembras como machos) de 4 y 18 meses de edad. Utilizaremos para el análisis datos de 5 hembras de 4 meses, 5 machos de 4 meses, 5 hembras de 18 meses y 5 machos de 18 meses, aunque para el ejemplo del pseudoalineamiento sólo usaremos un ratón.

## La estrategia

Para hacer el análisis lo más rápido posible, usaremos [salmon](https://combine-lab.github.io/salmon/). Salmon es un pseudoalineador que ha demostrado ser muy eficiente y comparable con los alineadores clásicos. Además, es rapidísimo. Posteriormente usaremos [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) para el análisis de expresión diferencial.

## 1: Pseudoalineamiento

### 1.1: Generar el index

Los indexes de salmon son muy fáciles de generar. Lo único que necesitamos son las secuencias en fasta de los **transcritos** a los que queremos (pseudo)alinear. Encontrarás dicho archivo en la carpeta */index* con el nombre *mm39-gencode-M36-transcript.fa.gz*. Estos datos se descargaron de gencode usando [el script gencode annotation](https://github.com/david-valle/gencode_annotation).

Para generar el index, corremos salmon en la terminal:

``` shell
salmon index -p 5 --gencode -t ./index/mm39-gencode-M36-transcript.fa.gz -i ./index/mm39.gencode.M36.salmon
```

Este comando usará 5 cores para generar el index y lo guardará en la carpeta *./index* con el nombre *mm39.gencode.M36.salmon*

Debe correr en menos de 10 minutos.

### 1.2: Hacer el pseudoalineamiento

Una vez creado el index, tenemos todo lo que necesitamos para realizar un pseudoalineamiento con salmon. Usamos el siguiente comando para (pseudo)alinear los archivos de la carpeta *fastq.* Recuerda que estos son sólo dos archivos que corresponden a una sola muestra. El dataset original tiene 20 muestras, así que para correrlo todo, el comando se va iterando sobre todos los archivos.

``` bash
salmon quant -p 5 -i ./index/mm39.gencode.M36.salmon -l A --seqBias --gcBias --discardOrphansQuasi -g ./index/mm39-gencode-M36-transcript_id-gene_id.txt -1 ./fastq/C57-cortex-4m-female-1_R1.fastq.gz -2 ./fastq/C57-cortex-4m-female-1_R2.fastq.gz -o mm39-C57-cortex-4m-female-1-salmon
```

El comando usará 5 cores para generar los pseudoalineamientos. Nota que usa una referencia adicional con la opción -g que es un archivo que contiene en una columna los ids de los transcritos y en otra los de los genes (y que se genera con el[script gencode annotation](https://github.com/david-valle/gencode_annotation)). Esto es para que el programa nos entregue una cuenta de genes, no sólo de transcritos. El resultado se guardará en la carpeta *mm39-C57-cortex-4m-female-1-salmon* y tomará aproximadamente 5 minutos o menos.

Dentro de los varios archivos que nos genera salmon, el que nos interesa se llama *quant.genes.sf*.

### 1.3: Generar la matriz de conteos

Para este paso, asumiremos que ya corrimos salmon en los 20 pares de archivos. Como no lo hemos hecho realmente, no corran los siguientes comandos, pero se los dejo de referencia. Esta serie de comandos nos permitirá generar una tabla de conteos en el folder counts. Ustedes ya tienen las tablas que resultan de correr estos comandos gracias a la magia de la televisión.

``` bash
cut -f 1,2 mm39-C57-cortex-4m-male-1-salmon/quant.genes.sf > counts/mm39-salmon-gene-annotation.tsv
cut -f1 mm39-C57-cortex-4m-male-1-salmon/quant.genes.sf > counts/mm39-C57-cortex-gene-count.tsv
for AGE in 4m 18m
do
    for SEX in male female 
    do
        for MICE in 1 2 3 4 5
        do
            cp mm39-C57-cortex-gene-count.tsv temp1
            cut -f5 mm39-C57-cortex-$AGE-$SEX-$MICE-salmon/quant.genes.sf | sed "s/NumReads/${SEX}_${AGE}_${MICE}/" > temp2
            paste temp1 temp2 > counts/mm39-C57-cortex-gene-count.tsv
            rm temp1 temp2
        done
    done
done
```

## 2: Análisis de expresión diferencial

### 2.1 Cargar los datos

Abran una sesión de Rstudio en la carpeta donde corrieron salmon. Lo primero que tenemos que hacer es cargar los datos de las counts a un objeto de tipo DGE_List que es un objeto específico de edgeR sobre el que haremos todo el análisis.

``` r
library(edgeR)

# Creamos los factores que vamos a ocupar en el análisis que nos dicen qué es cada cosa:
samples = factor(c(rep("m4",10),rep("m18",10)), levels=c("m4","m18"))
sex = factor(c(rep("male",5), rep("female",5),rep("male",5), rep("female",5)))
sample_col = c(rep("lightblue",5), rep("lightpink",5),rep("blue",5), rep("red",5))

# Cargamos las dos tablas en counts
gene_counts = read.delim("counts/mm39-C57-cortex-gene-count.tsv", row.names=1)
gene_annot = read.delim("counts/mm39-salmon-gene-annotation.tsv")

# Generamos nuestro objeto DGEList al que llamanos DGEList_age
DGEList_age = DGEList(counts=gene_counts, group=samples, genes=gene_annot)
```

### 2.2 Eliminar los genes con baja o nula expresión

Siempre es importante eliminar los genes poco expresados para que el análisis estadístico sea robusto. EdgeR tiene una función específica para ello. Adicionalmente, podemos ver con cuántos genes nos quedamos al final y guardar el dato en la variable n_features. Retuvimos un total de 17,315 genes.

``` r
keep = filterByExpr(DGEList_age)
n_features = sum(keep)
DGEList_age = DGEList_age[keep, keep.lib.sizes=FALSE]
```

### 2.3 Generar la matriz de diseño y calcular factores de normalización y estimar la dispersión.

La matriz de diseño le dice al programa qué clase de experimento tenemos. En este caso, simplemente cargamos las muestras. Posteriormente, calculamos los factores de normalización y dispersión que nos permiten realizar el análisis estadístico de expresión diferencial

``` r
design = model.matrix(~0+samples)
colnames(design) = c("m4","m18")

DGEList_age = calcNormFactors(DGEList_age)
DGEList_age = estimateDisp(DGEList_age, design=design)
```

### 2.4 Analizar el MDS plot y determinar el efecto de posibles confusores

Un MDS plot es muy similar a un PCA. Nos da esencialmente la misma información. Al inspeccionarlo notamos que los genes no sólo se separan por edad, sino también por sexo. EdgeR puede corregir por este factor confusor si lo consideramos relevante. Para analizar cómo quedarían los datos después de dicha corrección usamos la función removeBarchEffect.

``` r
# MDS plot con todos los datos
plotMDS(DGEList_age, col=sample_col)

# MDS plot eliminando el batch effect del sexo
plotMDS(removeBatchEffect(DGEList_age, batch=sex), col=sample_col)
```

Vemos que el tomar en cuenta el efecto del sexo no nos ayuda mucho, por lo que proseguiremos con el análisis sin realizar ninguna corrección adicional.

### 2.5 Preparar los contrastes y definir los cut-offs

Primero debemos determinar qué contrastes o comparaciones hacer. En este caso, dado que sólo tenemos una condición biológica relevante (la edad), haremos sólo una comparación. Pero en datasets más complejos se pueden hacer varias.

Adicionalmente, definimos nuestros cut-offs. En este caso usaremos un cut-off por FDR de 1e-5, y uno de log2 fold-chage mayor a 1.

Finalmente, hacemos el ajuste usando el método recomendado por edgeR (los modelos de quasi-likelihood GLM)

``` r
contrasts = makeContrasts(age=m18-m4, levels=design)
FDR = 1e-5
FC = 1
fit = glmQLFit(DGEList_age, design)
```

### 2.6 Prueba F para determinar genes up y down-regulated.

En el siguiente comando, determinamos cuáles son los genes que suben o que bajan según nuestros parámetros y guardamos los datos en tablas.

``` r
# Hacemos el F test para el contraste relevante
qlf <- glmQLFTest(fit, contrast=contrasts[,"age"])
# Guardamos todas las tags con la información de la comparación
DEG_qlf <- as.data.frame(topTags(qlf, n=n_features))

# Creamos filtros para guardar tablas sólo con los genes up- y down-regulated
up = (DEG_qlf$logFC > FC) & (DEG_qlf$FDR < FDR)
sum(up)
write.table(DEG_qlf[up,], paste("edgeR-DEG-age-Up-", FDR, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

down = (DEG_qlf$logFC < -FC) & (DEG_qlf$FDR < FDR)
sum(down)
write.table(DEG_qlf[down,], paste("edgeR-DEG-age-Down-", FDR, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
```

### 2.7 Volcano plot

Con el siguiente código, generamos un bonito volcano plot usando la función de plot de R. También se podría utilizar ggplot2 si les gusta más.

``` r
# Estos parámetros nos permiten modificar fácilmente los límites del volcano plot.
limx = 10
limy = 25

# Ploteamos los resultados separando los genes en aquellos que no cambian (grises), aquellos que suben (rojos) y aquellos que bajan (azules)
plot(DEG_qlf$logFC[!(up | down)], -log10(DEG_qlf$FDR[!(up | down)]), pch=19, col="gray", cex=0.4, xlab="log2 Expression fold change", ylab="-log10 FDR", main=paste("Volcano plot ", name, sep=""), xlim=c(-limx,limx),ylim=c(0,limy))
# Agregamos los up
points(DEG_qlf$logFC[up], -log10(DEG_qlf$FDR[up]), pch=19, col="red", cex=0.4)
# Agregamos los down
points(DEG_qlf$logFC[down], -log10(DEG_qlf$FDR[down]), pch=19, col="blue", cex=0.4)
# Agregamos las líneas que representan nuestros cutoffs
abline(h=-log10(FDR), col="black", lty=3)
abline(v=c(-FC,FC), col="black", lty=3)
# Guardamos todo
dev.copy(pdf,paste("Volcano_plot-RNA_seq-age-" FDR, ".pdf", sep=""))
dev.off()
```

### 2.8 Heatmap

Siempre es útil generar un heatmap para visualizar los resultados. Para ello, primero calcularemos la TPM a partir de la fpkm que edgeR genera. Esta es muy útil para otro tipo de plots también y nos puede servir para comparar el nivel de expresión incluso con otros datasets.

``` r
# Calculamos el fpkm con las funciones de edgeR. Noten que en realidad guardamos el log2(fpkm)
log2_fpkm_age = rpkm(DGEList_age, log=T)
log2_fpkm_age_average = rpkmByGroup(DGEList_age, log=T)

# Con esta función podemos pasar de FPKM a TMP
fpkm2tpm_log2 <- function(fpkm) { fpkm - log2(sum(2^fpkm)) + log2(1e6) }

# Aplicamos la función para generar los valores de TPM
log2_tpm_age = apply(log2_fpkm_age, 2, fpkm2tpm_log2)
log2_tpm_age_average = apply(log2_fpkm_age_average, 2, fpkm2tpm_log2)
```

Ahora, calculamos el z-score sólo para los genes que cambian significativamente y con esos datos generamos el heatmap

``` r
# Guardamos los ids de los genes significativos
id_significant = DEG_qlf$Name[up | down]
# Con esos ids hacemos un subset y calculamos el Z-score
zscore_significant = t(scale(t(log2_tpm_age[id_significant,]))) 

# Generamos un heatmap con el paquete complex heatmap
library(ComplexHeatmap)
Heatmap(zscore_significant, cluster_rows = T, cluster_columns = F, show_row_names = F, name = "Z score", km = 2, column_title = "All significant genes", row_title = c("Down","Up"))
```

Ese heatmap se ve muy bien, pero dado que son muchos genes, es imposible ponerles nombres. Hagamos un heatmap de sólo el top 20 de genes que cambian y además sustituyamos los ids por nombres de genes.

``` r
# Guardamos sólo el top 20 y calculamos el Z-score
DEG_top <- as.data.frame(topTags(qlf, n=20))
id_significant_top = DEG_top$Name
zscore_top = t(scale(t(log2_tpm_age[id_significant_top,]))) 

# Usamos un archivo que tiene los ids de los genes en una columna y los nombres en la otra para generar un diccionario que usaremos para poner nombres en lugar de ids. Este archivo se genra con el script de gencode annotation
gene_name_map = read.delim("./index/mm39-gencode-M36-gene_id-gene_name.txt", header=FALSE, row.names=1)

# Noten que ahora usamos como labels los nombres de los genes
Heatmap(zscore_top, cluster_rows = T, cluster_columns = F, row_labels = gene_name_map[id_significant_top,], name = "Z score", km = 2, column_title = "Top 20 significant genes", row_title = c("Down","Up"))
```

## 3. Análisis adicionales

Ya que tenemos la lista de genes que suben y que bajan, conviene saber un poco más de ellos. Para eso, siempre recomiendo hacer un análisis de Enriquecimiento con la plataforma [Enrichr](https://maayanlab.cloud/Enrichr/) y un análisis de interacción de proteínas usando [STRING](https://string-db.org/).

### 3.1 DAVID

Para que el análisis de GO de Enrichr sea robusto, necesitamos usar una lista de background relevante. En el caso de los estudios de RNA-seq el background debe ser la lista de genes que están expresados. Podemos obtener esa lista fácilmente en R, al igual que las ids de los genes up y downregulated:

``` r
# Enrichr sólo acepta nombres de genes. Así que obtenemos los nombres a partir de nuestro gene_name_map
id_all = DEG_qlf$Name
id_up = DEG_qlf$Name[up]
id_down = DEG_qlf$Name[down]

write.table(gene_name_map[id_all,], "expressed-genes-names.txt", quote=FALSE, row.names=F, col.names=F)
write.table(gene_name_map[id_up,], "upregulated-genes-names.txt", quote=FALSE, row.names=F, col.names=F)
write.table(gene_name_map[id_down,], "downregulated-genes-names.txt", quote=FALSE, row.names=F, col.names=F)
```

Noten que para tener una interpretación biológica más clara, es importante hacer el análisis separando a los genes up, de los down.

### 3.2 STRING

Para el análisis de STRING pueden jugar mucho con los tipos de evidencia y la confidence de las interacciones. Éstos van a tener un efecto significativo en las redes que generen. De igual manera, recomiendo hacer el análisis por separado, usando los up y down-regulated genes