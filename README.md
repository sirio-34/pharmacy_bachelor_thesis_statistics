# Bibliometric Analysis: Pharmaceutical Responses to CBRN Incidents

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Language](https://img.shields.io/badge/Language-R-276DC3.svg)](https://www.r-project.org/)

This repository contains the R script used for the publication analysis in my Pharmacy undergraduate thesis at Universidad Europea de Madrid.

The project examined pharmaceutical responses to five historical chemical, biological, radiological, and nuclear (CBRN) incidents. The publication analysis also included four comparison events.

The script retrieves yearly publication counts from PubMed through the NCBI Entrez API and generates statistical results, tables, and figures.

## My contribution and AI assistance

I defined the analysis requirements, including the PubMed searches and study periods, and reviewed the publication counts, tables, and figures. Gemini generated the R code based on my instructions.

## Methodology

The analysis compares publication activity before and after each selected event, excluding the event year. Most cases use eight-year periods before and after the incident; Salisbury uses seven-year periods.

Two statistical approaches are implemented:

1. **One-sided exact binomial test:** assesses post-event publication counts against a null proportion based on overall PubMed publication volume.
2. **Quasi-Poisson regression:** models yearly publication counts using a time term and the logarithm of total annual PubMed publications as an offset.

When the total pre-event count is zero, the script adds `0.25` to the annual incident-related counts. It also applies a minimum dispersion value of `1` when calculating standard errors.

These adjustments do not establish causality or guarantee protection against false positives.

## Requirements

The script requires R and an internet connection. It automatically attempts to install missing packages:

- `rentrez`: access to NCBI Entrez.
- `ggplot2`, `dplyr`, `tidyr`, and `scales`: data handling and visualization.

## NCBI configuration

An NCBI API key is optional. The default request limits are three requests per second without a key and ten with a key. See the [NCBI usage guidelines](https://www.ncbi.nlm.nih.gov/books/NBK25497/pdf/Bookshelf_NBK25497.pdf).

The script contains a placeholder in `set_entrez_key()`. Replace it locally with your own key:

```r
rentrez::set_entrez_key("YOUR_NCBI_API_KEY")
```

If running without a key, comment out or remove the placeholder call.

The `pausa` parameter in `consultar_pubmed()` controls the delay between requests. Its current default is `0.6` seconds. Keep requests within NCBI's usage limits.

## Interpretation

This is an exploratory analysis of publication activity. Results depend on the search terms, study periods, and PubMed records retrieved. Review query warnings and missing values before interpreting the outputs.

Changes in publication counts do not establish that an incident caused a change in research activity.

## License

Released under the MIT License. See `LICENSE`.

## Original Spanish documentation

The Spanish text below preserves the README used when this repository was assessed as part of my Pharmacy undergraduate thesis. The English section above was added later for international readers and includes clarifications about the methodology and configuration.

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
```
