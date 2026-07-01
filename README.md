# Análisis Bibliométrico: Respuesta Farmacológica a Incidentes NRBQ

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Language](https://img.shields.io/badge/Language-R-276DC3.svg)](https://www.r-project.org/)

Repositorio oficial con el código fuente utilizado para el Trabajo de Fin de Grado (TFG) sobre el impacto bibliométrico y la respuesta farmacológica y toxicológica ante incidentes NRBQ (Nuclear, Radiológico, Biológico y Químico). 

Este script en R automatiza la extracción de datos desde la base de datos médica PubMed (vía la API de NCBI Entrez) y ejecuta un modelado estadístico riguroso para medir el impacto real de desastres históricos en la literatura científica.

## Metodología

El diseño de investigación utiliza un enfoque de **Casos y Controles Históricos**. Compara el volumen de publicaciones indexadas 8 años antes (Pre) frente a 8 años después (Post) del incidente. 

Para aislar el impacto del desastre del crecimiento natural (vegetativo) de la propia base de datos, el script ejecuta un análisis dual:
1. **Test Binomial Exacto (unilateral):** Para contrastar el incremento absoluto frente a la proporción nula esperada.
2. **Regresión de Quasipoisson (con ajuste de Offset):** Un modelo lineal generalizado que utiliza el logaritmo de las publicaciones totales anuales de PubMed como `offset`. Incluye una corrección de continuidad de `+0.25` para evitar errores de separación perfecta en casos de censura o ausencia previa de literatura (ej. el agente Novichok antes de 2018), y aplica un suelo de dispersión para blindar a los controles frente a falsos positivos.

## Requisitos y Configuración

El script está diseñado para instalar automáticamente las dependencias si no se encuentran en el sistema:
* `rentrez` (Conexión a la API de NCBI)
* `ggplot2`, `dplyr`, `tidyr`, `scales` (Manipulación de datos y visualización)

**Importante:** Para evitar bloqueos por parte de los servidores gubernamentales de EE.UU., es necesario introducir tu propia clave API de NCBI en la línea 24 del código, si no tienes que ponerle mas tiempo de espera de entre peticiones en la línea 84, con una pausa de 0.8 iria bien:
```R
set_entrez_key("xxxxxxxxxxxxxxxx")
consultar_pubmed <- function(query, year, pausa = X.X)
