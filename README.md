# Alineamiento de secuencias de ADN
Se suministra a los alumnos un código secuencial (ver apéndices al final de la página) que utiliza un método de fuerza bruta para detectar alineamientos (coincidencias) exactas de una serie de patrones en una cadena genética más larga.

Ejemplo:
Secuencia original: CCGTACCTGTGACCGTTAAAACTTTC
Patrón: ACCGT
Alinea en posición 11 (la primera posición es la 0)

El programa simula una situación real generando según los argumentos suministrados una secuencia aleatoria, y una serie de patrones que pueden ser o completamente aleatorios (puede que se encuentren o no) o muestras escogidas aleatoriamente de la secuencia original (siempre se acabarán encontrando).

El programa recorre la lista de patrones buscando una coincidencia exacta a partir de cada una de las posibles posiciones de comienzo en la secuencia original. Para la búsqueda en la primera coincidencia. Hay algoritmos más eficientes, pero vamos a utilizar exclusivamente este método de fuerza bruta porque es más sencillo, muy regular y paralelizable, facilitando probar los conceptos de la asignatura.

Resultados

El programa calcula además la cantidad total de patrones que se encuentran, para cada posición de la cadena original la cantidad de patrones que alinean sobre esa posición, y el patrón de longitud más larga que alinea en esa posición. Realiza y muestra en la salida por defecto una serie de checksums para verificar fácilmente el resultado comparando con los resultados del programa secuencial.

# Detalles del código secuencial
Argumentos del programa

    * <seq_length>: Longitud de la secuencia proincipal
    * <prob_G>: Probabilidad de aparición de nucleótidos G
    * <prob_C>: Probabilidad de aparición de nucleótidos C
    * <prob_A>: Probabilidad de aparición de nucleótidos A
    * (La probabilidad de aparición de nucleóticos T es uno menos la suma de las tres anteriores)
    * <pat_rng_num>: Número de patrones aleatórios
    * <pat_rng_length_mean>: Longitud media de esos patrones
    * <pat_rng_length_dev>: Desviación de la longitud de esos patrones
    * <pat_samples_num>: Número de patrones que son muestras de la secuencia original
    * <pat_samp_length_mean>: Longitud media de esos patrones 
    * <pat_samp_length_dev>: Desviación de la longitud de esos patrones
    * <pat_samp_loc_mean>: Localización media del comienzo de las muestras
    * <pat_samp_loc_dev>: Desviación de la localización de comienzo
    * <pat_samp_mix:B[efore]|A[fter]|M[ixed]>: Esta parámetro indica:
    * B: Los patrones de muestras están en la lista delante de los aleatorios
    * A: Los patrones de muestras están en la lista detrás de los aleatorios
    * M: Los patrones de muestras y aleatorios están entremezclados
    * <long_seed>: Una semilla aleatoria para la generación de casos diferentes y reproducibles para los mismos argumentos.
    * Estos argumentos dan al usuario mucho control para la generación de escenarios de diferentes tipos y la colocación de la carga computacional. Se puede utilizar el modo DEBUG con tamaños pequeños para hacerse una idea de lo que se puede obtener con diferentes combinaciones.

## Objetivo

El objetivo es paralelizar y optimizar el código sin modificar el algoritmo que se aplica, ni los resultados de la simulación. Es importante distinguir los trozos de código que pueden ser paralelizados directamente, los que necesitan modificaciones, resolver condiciones de carrera, etc.

## Resultados

El programa muestra el tiempo de ejecución de la parte de computo (sin la inicialización), y los checksums que se utilizan para verificar la corrección de la simulación.

## Paralelizacion y Optimizacion
Parte del programa que se ha paralelizado es la indicada  "START HERE:"  Claramente indicado en los comentarios de cada programa
La idea del algoritmo es resolver mediante el uso de fuerza bruta. Ese concepto no es alterado en ninguna de las versiones paralelizadas.

Tipos de paralelizacion usada:
    * Sequencial (align_seq) => Version secuencial del algoritmo de fuerza bruta
    * OpenMP (Open MultiProcessing) API (align_omp.c) => Simple paralelizacion de bucles y creacion the threads para resolver las tareas mas pesadas
    * MPI (Message Passing Interface) standard (align_mpi) => Uso de diferentes procesos en diferentes maquinas para paralelizar tareas pesadas computacionalmente mediante el paso de mensajes entre procesos. Resolucion mediante patron Master-Slave un procesos se encarga de mandar a los demas las tareas y unifica los resultados.
    * CUDA (Computed Unified Device Architecture) API/NVIDIA => Paralelizacion mediante el uso GPU aprovechando el gran numero de cores disponibles para abordar la tarea de fuerza bruta de una forma paralela mediante el uso de funciones kernel y uso memoria dinamica en la CPU

Uso de librerias: CUDA/NVIDIA, MPICH, OpenMP.

# Pruebas 

Las pruebas iniciales se realizaron sobre el HPC (High Performance Computing) cluster del grupo de investigacion Trago de la Universidad de Valladolid del Grado de Ingenieria Informatica.

- Maquina Gorgon (OpenMP, MPI) => 2 x AMD EPYC 7713 64-Core Processor @ 2Ghz (total 128 cores / 256 hilos) & 1 NVIDIA A100 80GB, 13824 cores @ 1.4Mhz.
- Maquina Ampere (CUDA) => 1 NVIDIA A100 80GB, 13824 cores @ 1.4Mhz

Ejemplo de ejecución en local (nvidia y mpi necesitaran sus entornos especificos instalados aparte de la librerias):
```
./align_type 500000 0.15 0.4 0.15 1120 20 10 2110 50000 500 7000 4000 M 346536
```
# Compilación y máquinas de referencia
Para medir la eficiencia el sistema tablón y sus leaderboard compilan los códigos con un compilador GCC (v10.3), con -O3 como único flag de optimización automática. La ejecución se realiza en máquinas de referencia que se detallan en el FAQ de tablón para cada práctica.

```
$ make all => // To generate all binaries OpenMP, MPI, CUDA, Sequencial
$ make align_type => // To generate only one type of binary
$ make clean => // Remove all the binaries and object files
```
## Modo Debug

Compilar el el programa con -DDEBUG (lo puede hacer el makefile suministrado automáticamente), se activan unas partes del programa que muestran los argumentos leidos, las cadenas generadas y los resultados de los arrays auxiliares que se utilizan para los checksums.

