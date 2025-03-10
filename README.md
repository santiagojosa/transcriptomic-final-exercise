## Ejercicio final de transcriptómica (curso 2024-2025)

El objetivo de este trabajo es que el alumno demuestre los conocimientos obtenidos acerca del análisis de RNA-seq. Para ello deberá **redactar un informe** en el que se expliquen los datos de partida y se extraiga una **conclusión de los resultados**. Así mismo, el alumno debe **detallar el proceso de análisis** indicando el *software* (incluída la versión) empleado, así como los parámetros utilizados en cada uno de los pasos. En caso de que hubiera que eliminar muestras por motivos técnicos o biológicos, el alumno debe indicar y justificar el por qué en cada caso. Se pueden introducir *code chunks*, imágenes  y tablas para apoyar el informe. De manera alternativa, puede aportarse todo el código en forma de repositorio público. Se han planteado **5 preguntas (10 puntos en total)** para guiar la redacción del informe, que debe tener una **extensión máxima de 8.000 palabras, sin contar pies de figura o tablas.**

El trabajo consta de dos apartados en los que se utilizarán datos de un experimento de 24 cultivos primarios de tumores paratiroideos negativos para receptores de estrógenos alfa (ERα). Las muestras, procedentes de 4 pacientes diferentes, se han tratado con dos fármacos diferentes: diarilpropionitrilo (DPN) o 4-hidroxitamoxifeno (OHT) a 24h o 48h. El DPN es un agonista del ERα mientras que el OHT es un inhibidor competitivo de los receptores de estrógenos.

El **primer apartado (3 preguntas)** abarca los pasos de control de calidad y de fuentes de contaminación, *trimming*, alineamiento y cuantificación para obtener cuentas crudas y normalizadas a partir de un **subset de ficheros fastq**. El **segundo apartado (2 preguntas)** parte de la **matriz completa de cuentas crudas** y está enfocado en realizar un control de calidad biológico, detectar los genes diferencialmente expresados entre condiciones y los pathways enriquecidos en cada una de ellas.

### Apartado 1
El dataset original consta de 27 muestras paired-end depositadas en SRA. Con el fin de poder abordar las cuestiones planteadas a en el primer apartado sólo es necesario descargar 2 muestras (SRR479052 y SRR479054), es decir cuatro ficheros fastq. Además, en este repositorio se proporciona:

- Un fichero fasta con la secuencia de la referencia genómica, en este caso correspondiente al cromosoma 21 humano (ensamblaje GRCh38).
- Un fichero GTF con la anotación génica para los genes del cromosoma 21 (GRCh38.ensembl.109).

Se pide realizar un análisis de dichas muestras similar al realizado en clase, considerando los siguientes puntos:

**Pregunta 1 (1.5 puntos):** Realizar control de calidad de dichas muestras con el programa FastQC, incluir plots más reseñables y comentar cada uno de los apartados. De manera complementaria se podrá realizar un analisis de contaminación con el programa FastQScreen.

**Pregunta 2 (1.5 puntos):** Para poder llevar a cabo el alineamiento de las muestras en vuestros ordenadores será necesario trabajar con archivos reducidos correspondientes al cromosoma 21. Se requiere el indexado de la secuencia de este cromosoma, así como el alineamiento de las muestras a dicha referencia. Para ello se podrá utilizar el alineador HISAT2 utilizado en clase u otros (alineadores o pseudoalineadores). Comentar cada uno de los comandos y parámetros empleados, justificando su uso.

**Pregunta 3 (1.5 puntos):** Una vez generados los archivos alineados se reportarán las estadísticas de alineamiento y se procederá a la cuantificación de la expresión utilizando el archivo GTF correspondiente. Para ello se podrá utilizar HTSeq u otras herramientas de cuantificación. En cualquier caso, detallar y justificar los comandos y parámetros empleados para ello.

### Apartado 2
En este repositorio se proporcionan todos los inputs del Apartado 2:

- La matriz de cuentas crudas de los 24 cultivos analizados. 
- Data frame con los metadatos asociados al experimento.
- GMT para realizar un GSEA.

**Pregunta 4 (3 puntos):** ¿Qué genes se encuentran diferencialmente expresados entre las muestras pertenecientes al grupo tratado con OHT con respecto al control tras 24h? ¿Y en el caso de las muestras tratadas con DPN, también tras 24h? Como parte de la respuesta a esta pregunta, podéis entregar una o más tablas adjuntas donde se incluyan los genes diferencialmente expresados, indicando el criterio o los criterios que habéis seguido para filtrar los resultados, así como cualquier otro gráfico o gráficos que hayáis elaborado durante el análisis.

**Pregunta 5 (2.5 puntos):** Nuestro colaborador ha comparado las muestras tratadas con DPN tras 48h con las muestras control. Nos ha llamado para contarnos que los cambios de expresión tras 48h son mucho más evidentes y nos preguntamos si el DPN produce algún efecto en las primeras 24h. Para contestar esta pregunta, le pedimos que genere un GMT (disponible en la carpera de input) con los genes más expresados en las muestras tratadas tras 48h (DPN_perturbed) y los genes más expresados en la muestra control (DPN_unperturbed). Realiza un análisis con GSEA para determinar el efecto del tratamiento a las 24h. ¿A qué conclusión llegas? Incluye una tabla con los resultados del análisis, destacando las columnas con los valores utilizados para extraer vuestras conclusiones. También incluye los gráficos característicos de este tipo de análisis.
