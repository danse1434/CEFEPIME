# Farmacocinética Poblacional de Cefepime 

Este repositorio se crea con el fin de preservar las funciones del código R para el mantenimiento del mismo, así como mejoras en legibilidad y rendimiento. 

## Introducción
> La presencia de periodos de neutropenia prolongados es frecuente en pacientes que están recibiendo tratamiento quimioterapéutico para el cáncer. La presencia de estos episodios de neutropenia, dejan al paciente con un riesgo > importante de presentar infecciones severas, que pueden tener en morbilidad o mortalidad en estos pacientes. Cuando se presenta fiebre durante estos periodos de neutropenia, se debe pensar como un síntoma infección a menos > que se pueda demostrar lo contrario. 

Al presentarse una infección se deben realizar procedimientos diagnósticos tendientes a la identificación del patógeno causante, así como una evaluación física completa que permita identificar algún foco de infección. Esto en conjunto con la administración de una terapia antimicrobiana que tenga la mayor efectividad posible dentro de un tiempo prudente. La efectividad de la terapia se puede estimar mediante cálculos a partir de la farmacocinética y farmacodinamia del antimicrobiano. Se puede establecer la probabilidad de alcanzar un objetivo PK-PD y establecer lineamientos de dosis incorporando factores del microorganismo como la CMI. 

La variabilidad en la respuesta a un tratamiento depende de varios factores, pero posiblemente uno de los que más interfiere - en estos pacientes – corresponde a los cambios fisiopatológicos a los que se ven sometidos estos pacientes. Los parámetros farmacocinéticos en los pacientes con NF pueden ser alterados por una amplia variedad de condiciones y esta alteración puede ser la causa de una menor efectividad en los tratamientos. Pese a la magnitud del problema, no se conoce de manera integral el comportamiento de muchos antimicrobianos en personas con NF posquimioterapia. Los enfoques que implican la integración de modelos de enfermedad, con parámetros farmacocinéticos y farmacodinámicos han ganado amplia aceptación en tiempos recientes, ya que estos pueden optimizar la respuesta terapéutica.

Un objetivo de este proyecto de investigación es determinar los parámetros farmacocinéticos de cefepime, piperacilina y vancomicina en pacientes con neoplasias hematológicas y neutropenia febril posquimioterapia. Lo anterior, tiene como fin optimizar la terapia que se les brinda a estos pacientes. El presente estudio se trataría de un análisis post-hoc de los datos obtenidos como parte de una investigación en esta población atendida en el Instituto Nacional de Cancerología (INC). 

Los parámetros PK serán obtenidos mediante la construcción de un modelo basado en farmacocinética poblacional, con esto se podrá conocer el valor medio y variabilidad de los parámetros. El comportamiento de los mismos frente a factores fisiopatológicos relevantes del paciente; la proporción de variabilidad observada debida a diferencias entre individuos; entre otros. 

Los resultados provenientes de este estudio podrían ser aplicados en clínica con el fin de conocer la dosis en pacientes individuales. Entre las limitaciones se reconoce que este proyecto no contempla la realización de una validación externa, sin embargo, se realizarán procedimientos de evaluación interna de referencia en la literatura. 

El estudio sentará un precedente para la futura realización de estudios de farmacocinética poblacional en pacientes colombianos. Este enfoque de análisis es más eficiente en el modelamiento de datos farmacocinéticos, pues permite obtener información que normalmente es ignorada. Lo anterior con integración entre covariables propias del individuo y la disposición de los fármacos en el organismo. Este tipo de modelos podrían ser usados en la práctica clínica, en virtud de que se han desarrollado metodologías que permiten integrarlos con monitorización de fármacos. Al obtener un punto de concentración en el tiempo en un paciente, se pueden estimar los parámetros farmacocinéticos del individuo mediante análisis bayesiano y relacionarlos con la probabilidad de un resultado favorable o desfavorable.

## Contenido
Capítulo | Contenido
------------ | -------------
00 | Manipulación inicial de los datos (cierre). Reporte de validación de metodología analítica. 
01 | Análisis Exploratorio de los Datos (EDA)
02 | Desarrollo de modelo farmacocinético base
03 | Exploración del efecto de la censura de los datos
04 | Determinación de modelo de error residual. 
05 | Remuestreo por bootstrap para el modelobase
06 | Evaluación de correlación entre etas.
07 | Evaluación de propiedades de convergencia del modelo. 
08 | Verificación de outliers y análisis de influencia.
09 | Paso COV
10 | Refinamiento del modelo base
11 | Evaluación exploratoria de las covariables.
12 | Evaluación de modelo de covariables
13 | Selección de modelo final de covariables 1
14 | Selección de modelo final de covariables 2
15 | Descripción del modelo final de covariables
97 | Comparativos entre modelos

