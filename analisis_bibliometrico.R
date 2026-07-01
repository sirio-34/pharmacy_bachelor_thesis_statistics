###############################################################################
#  ANÁLISIS BIBLIOMÉTRICO – TFG Respuesta Farmacológica a Incidentes NRBQ
#  -----------------------------------------------------------------------
#  Autor : Diego
#  Fecha : Junio 2026
#  Descripción:
#    Compara el volumen de publicaciones en PubMed 8 años antes vs. 8 años
#    después de incidentes NRBQ específicos, controlando por el crecimiento
#    vegetativo de PubMed. Aplica un Test Binomial Exacto unilateral y una
#    Regresión Quasipoisson con offset del total de PubMed.
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 0. INSTALACIÓN AUTOMÁTICA DE PAQUETES
# ─────────────────────────────────────────────────────────────────────────────
paquetes_necesarios <- c("rentrez", "ggplot2", "dplyr", "tidyr", "scales")

for (pkg in paquetes_necesarios) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0(">>> Instalando paquete: ", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(rentrez)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

set_entrez_key("here you put your API key from NCBI")

message("\n========================================")
message("   ANÁLISIS BIBLIOMÉTRICO NRBQ – TFG")
message("========================================\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. DEFINICIÓN DE CASOS Y CONTROLES
# ─────────────────────────────────────────────────────────────────────────────
casos <- list(
  # ── CASOS PRINCIPALES ──
  Bhopal = list(
    t_year = 1984,
    query  = "\"methyl isocyanate\"[tiab] AND (toxicity[tiab] OR poisoning[tiab] OR antidote[tiab] OR treatment[tiab] OR pulmonary[tiab] OR exposure[tiab])"
  ),
  Chernobyl = list(
    t_year = 1986,
    query  = "(\"radiation emergency\"[tiab] OR \"radiation protection\"[tiab] OR \"nuclear accident\"[tiab] OR \"nuclear accident\"[Mesh] OR \"thyroid blockade\"[tiab] OR prophylaxis[tiab]) AND (\"potassium iodide\"[Mesh] OR \"potassium iodide\"[tiab])"
  ),
  Goiania = list(
    t_year = 1987,
    query  = "(\"cesium\"[tiab] OR \"cesium-137\"[tiab]) AND (\"acute radiation syndrome\"[Mesh] OR \"radiation injuries\"[Mesh] OR treatment[tiab] OR \"prussian blue\"[tiab] OR radiogardase[tiab] OR decorporation[tiab])"
  ),
  Sverdlovsk = list(
    t_year = 1979,
    query  = "anthrax[Mesh] AND (vaccine[tiab] OR treatment[tiab] OR prophylaxis[tiab] OR antidote[tiab])"
  ),
  Salisbury = list(
    t_year = 2018,
    query  = "(\"Novichok\"[tiab] OR \"nerve agent\"[tiab] OR \"nerve agents\"[tiab]) AND (oxime[tiab] OR antidote[tiab] OR reactivation[tiab] OR treatment[tiab] OR atropine[tiab] OR therapy[tiab] OR poisoning[tiab])"
  ),
  # ── CONTRAPARTES HISTÓRICAS ──
  Anthrax_Contra = list(
    t_year = 2001,
    query  = "anthrax[Mesh] AND (vaccine[tiab] OR treatment[tiab] OR prophylaxis[tiab] OR antidote[tiab])"
  ),
  Fukushima_Contra = list(
    t_year = 2011,
    query  = "(\"radiation emergency\"[tiab] OR \"radiation protection\"[tiab] OR \"nuclear accident\"[tiab] OR \"nuclear accident\"[Mesh] OR \"thyroid blockade\"[tiab] OR prophylaxis[tiab]) AND (\"potassium iodide\"[Mesh] OR \"potassium iodide\"[tiab])"
  ),
  West_Virginia_Contra = list(
    t_year = 2008,
    query  = "\"methyl isocyanate\"[tiab] AND (toxicity[tiab] OR poisoning[tiab] OR antidote[tiab] OR treatment[tiab] OR pulmonary[tiab] OR exposure[tiab])"
  ),
  Mayapuri_Contra = list(
    t_year = 2010,
    query  = "(\"cobalt\"[tiab] OR \"cobalt-60\"[tiab]) AND (\"acute radiation syndrome\"[Mesh] OR \"radiation injuries\"[Mesh] OR treatment[tiab])"
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# 2. FUNCIÓN AUXILIAR: Consulta PubMed por año
# ─────────────────────────────────────────────────────────────────────────────
consultar_pubmed <- function(query, year, pausa = 0.6) {
  # Publicaciones específicas del tema
  query_year <- paste0("(", query, ") AND (\"", year, "\"[PDAT])")
  n_especifico <- tryCatch({
    res <- entrez_search(db = "pubmed", term = query_year, retmax = 0)
    res$count
  }, error = function(e) {
    message(paste0("  [AVISO] Error en consulta específica para ", year, ": ", e$message))
    Sys.sleep(2)
    # Reintento
    tryCatch({
      res <- entrez_search(db = "pubmed", term = query_year, retmax = 0)
      res$count
    }, error = function(e2) {
      message(paste0("  [ERROR] Fallo definitivo para ", year))
      return(NA)
    })
  })
  
  Sys.sleep(pausa)
  
  # Total de publicaciones en PubMed ese año
  query_total <- paste0("\"", year, "\"[PDAT]")
  n_total <- tryCatch({
    res <- entrez_search(db = "pubmed", term = query_total, retmax = 0)
    res$count
  }, error = function(e) {
    message(paste0("  [AVISO] Error en consulta total para ", year, ": ", e$message))
    Sys.sleep(2)
    tryCatch({
      res <- entrez_search(db = "pubmed", term = query_total, retmax = 0)
      res$count
    }, error = function(e2) {
      message(paste0("  [ERROR] Fallo definitivo total para ", year))
      return(NA)
    })
  })
  
  Sys.sleep(pausa)
  
  return(list(n_especifico = n_especifico, n_total = n_total))
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. RECOLECCIÓN DE DATOS AÑO POR AÑO
# ─────────────────────────────────────────────────────────────────────────────
datos_todos <- data.frame()

for (nombre in names(casos)) {
  caso   <- casos[[nombre]]
  t_year <- caso$t_year
  query  <- caso$query
  
  # Excepción de 7 años para Salisbury (evita el año incompleto 2026)
  if (nombre == "Salisbury") {
    anios_pre  <- (t_year - 7):(t_year - 1)
    anios_post <- (t_year + 1):(t_year + 7)
  } else {
    anios_pre  <- (t_year - 8):(t_year - 1)
    anios_post <- (t_year + 1):(t_year + 8)
  }
  
  todos_anios <- c(anios_pre, anios_post)
  
  message(paste0("────────────────────────────────────────"))
  message(paste0("▶ Procesando: ", nombre, "  (año incidente: ", t_year, ")"))
  message(paste0("  Pre:  ", min(anios_pre), "-", max(anios_pre)))
  message(paste0("  Post: ", min(anios_post), "-", max(anios_post)))
  message(paste0("────────────────────────────────────────"))
  
  df_caso <- data.frame(
    caso       = character(0),
    year       = integer(0),
    periodo    = character(0),
    n_pub      = integer(0),
    n_total    = integer(0),
    stringsAsFactors = FALSE
  )
  
  for (anio in todos_anios) {
    periodo <- ifelse(anio < t_year, "Pre", "Post")
    message(paste0("  Consultando ", anio, " (", periodo, ")..."))
    
    res <- consultar_pubmed(query, anio)
    
    df_caso <- rbind(df_caso, data.frame(
      caso    = nombre,
      year    = anio,
      periodo = periodo,
      n_pub   = as.integer(res$n_especifico),
      n_total = as.integer(res$n_total),
      stringsAsFactors = FALSE
    ))
  }
  
  datos_todos <- rbind(datos_todos, df_caso)
  message(paste0("  ✓ ", nombre, " completado.\n"))
}

message("\n✓ Recolección de datos finalizada.\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4. ANÁLISIS ESTADÍSTICO
# ─────────────────────────────────────────────────────────────────────────────
resultados <- data.frame()

for (nombre in names(casos)) {
  message(paste0("\n══════════════════════════════════════"))
  message(paste0("  RESULTADOS: ", nombre))
  message(paste0("══════════════════════════════════════"))
  
  df <- datos_todos %>% filter(caso == nombre)
  
  df_pre  <- df %>% filter(periodo == "Pre")
  df_post <- df %>% filter(periodo == "Post")
  
  sum_pub_pre   <- sum(df_pre$n_pub,   na.rm = TRUE)
  sum_pub_post  <- sum(df_post$n_pub,  na.rm = TRUE)
  sum_tot_pre   <- sum(df_pre$n_total, na.rm = TRUE)
  sum_tot_post  <- sum(df_post$n_total, na.rm = TRUE)
  
  # ── 4a. Test Binomial Exacto (unilateral: post > pre) ──
  total_pub     <- sum_pub_pre + sum_pub_post
  p0            <- sum_tot_post / (sum_tot_pre + sum_tot_post) # Proporción nula
  
  if (total_pub > 0 & !is.na(p0) & p0 > 0 & p0 < 1) {
    binom_test <- binom.test(
      x           = sum_pub_post,
      n           = total_pub,
      p           = p0,
      alternative = "greater"
    )
    binom_p   <- binom_test$p.value
    # IRR binomial: (post/total_post) / (pre/total_pre)
    rate_pre  <- sum_pub_pre  / sum_tot_pre
    rate_post <- sum_pub_post / sum_tot_post
    irr_binom <- ifelse(rate_pre > 0, rate_post / rate_pre, Inf)
  } else {
    binom_p   <- NA
    irr_binom <- NA
  }
  
  message(paste0("\n  [Test Binomial Exacto]"))
  message(paste0("    Pub Pre:  ", sum_pub_pre, "  |  Pub Post: ", sum_pub_post))
  message(paste0("    Total PubMed Pre:  ", sum_tot_pre))
  message(paste0("    Total PubMed Post: ", sum_tot_post))
  message(paste0("    Proporción nula (p0): ", round(p0, 4)))
  message(paste0("    IRR: ", ifelse(is.finite(irr_binom), round(irr_binom, 3), "Inf")))
  message(paste0("    p-valor: ", format(binom_p, digits = 4, scientific = TRUE)))
  message(paste0("    Significativo (<0.05): ", ifelse(!is.na(binom_p) & binom_p < 0.05, "SÍ ***", "NO")))
  
  # ── 4b. Regresión Quasipoisson ──
  df_model <- df %>%
    mutate(
      periodo_factor = factor(periodo, levels = c("Pre", "Post")),
      year_index     = year - min(year)
    )
  
  # Corrección de separación perfecta si Pre == 0
  aplicada_correccion <- FALSE
  if (sum_pub_pre == 0) {
    df_model$n_pub <- df_model$n_pub + 0.25
    aplicada_correccion <- TRUE
    message("    [Corrección aplicada: +0.25 a n_pub por Pre = 0]")
  }
  
  qpois_irr  <- NA
  qpois_p    <- NA
  qpois_disp <- NA
  
  tryCatch({
    modelo <- glm(
      n_pub ~ periodo_factor + year_index + offset(log(n_total)),
      data   = df_model,
      family = quasipoisson(link = "log")
    )
    
    coefs <- summary(modelo)$coefficients
    disp  <- summary(modelo)$dispersion
    # Suelo de dispersion: no deflactar por debajo de Poisson. Evita falsos
    # positivos cuando hay INFRADISPERSION (disp < 1) en series escasas, que es
    # lo que hacia "significativos" a controles como Fukushima o West Virginia.
    disp_floor <- max(1, disp)
    if ("periodo_factorPost" %in% rownames(coefs)) {
      beta <- coefs["periodo_factorPost", "Estimate"]
      se   <- coefs["periodo_factorPost", "Std. Error"] * sqrt(disp_floor / disp)
      tval <- beta / se
      qpois_irr <- exp(beta)
      # Una cola SOLO en los 5 casos reactivos (hipotesis direccional post > pre).
      # Dos colas en los controles, donde se contrasta la AUSENCIA de cambio
      # (no hay hipotesis direccional). Dividir p/2 en los controles inflaba su
      # significacion de forma indebida.
      es_contrap <- grepl("Contra", nombre)
      if (es_contrap) {
        qpois_p <- 2 * pt(-abs(tval), df = modelo$df.residual)
      } else {
        qpois_p <- pt(-tval, df = modelo$df.residual)
      }
    }
    qpois_disp <- disp
  }, error = function(e) {
    message(paste0("    [AVISO] Error en modelo Quasipoisson: ", e$message))
  })
  
  message(paste0("\n  [Regresión Quasipoisson]"))
  message(paste0("    IRR:       ", ifelse(!is.na(qpois_irr), round(qpois_irr, 3), "N/A")))
  message(paste0("    p-valor:   ", ifelse(!is.na(qpois_p), format(qpois_p, digits = 4, scientific = TRUE), "N/A")))
  message(paste0("    Dispersión:", ifelse(!is.na(qpois_disp), round(qpois_disp, 3), "N/A")))
  message(paste0("    Significativo (<0.05): ",
                 ifelse(!is.na(qpois_p) & qpois_p < 0.05, "SÍ ***", "NO")))
  if (aplicada_correccion) {
    message(paste0("    (Se aplicó corrección +0.25 por separación perfecta)"))
  }
  
  # Guardar resultados
  resultados <- rbind(resultados, data.frame(
    Caso                     = nombre,
    Anio_Incidente           = casos[[nombre]]$t_year,
    Pub_Pre                  = sum_pub_pre,
    Pub_Post                 = sum_pub_post,
    Total_PubMed_Pre         = sum_tot_pre,
    Total_PubMed_Post        = sum_tot_post,
    Binom_IRR                = irr_binom,
    Binom_pvalor             = binom_p,
    Binom_Significativo      = ifelse(!is.na(binom_p) & binom_p < 0.05, "SÍ", "NO"),
    QuasiP_IRR               = qpois_irr,
    QuasiP_pvalor            = qpois_p,
    QuasiP_Dispersion        = qpois_disp,
    QuasiP_Significativo     = ifelse(!is.na(qpois_p) & qpois_p < 0.05, "SÍ", "NO"),
    Correccion_Separacion    = aplicada_correccion,
    stringsAsFactors         = FALSE
  ))
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. EXPORTAR CSV
# ─────────────────────────────────────────────────────────────────────────────
archivo_csv <- "resultados_bibliometria_dual_TFG.csv"
write.csv(resultados, file = archivo_csv, row.names = FALSE, fileEncoding = "UTF-8")
message(paste0("\n✓ Resultados guardados en: ", archivo_csv))

# También guardamos los datos brutos por año
write.csv(datos_todos, file = "datos_brutos_por_anio_TFG.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
message("✓ Datos brutos por año guardados en: datos_brutos_por_anio_TFG.csv")

# ─────────────────────────────────────────────────────────────────────────────
# 6. GRÁFICO DE BARRAS AGRUPADAS
# ─────────────────────────────────────────────────────────────────────────────

# Orden: cada caso principal (por año, de más antiguo a más reciente) seguido
# de su contraparte; Salisbury al final por no tener contraparte.
orden_casos <- c(
  "Sverdlovsk", "Anthrax_Contra",        # 1979 / ántrax EE.UU.
  "Bhopal",     "West_Virginia_Contra",  # 1984 / West Virginia
  "Chernobyl",  "Fukushima_Contra",      # 1986 / Fukushima
  "Goiania",    "Mayapuri_Contra",       # 1987 / Mayapuri
  "Salisbury"                            # 2018 / sin contraparte
)

# Etiqueta de par (cada caso principal con su contraparte) y orden cronológico
niveles_par <- c(
  "Sverdlovsk / Ántrax EE.UU.",
  "Bhopal / West Virginia",
  "Chernóbil / Fukushima",
  "Goiânia / Mayapuri",
  "Salisbury"
)
grupo_par <- c(
  "Sverdlovsk"           = "Sverdlovsk / Ántrax EE.UU.",
  "Anthrax_Contra"       = "Sverdlovsk / Ántrax EE.UU.",
  "Bhopal"               = "Bhopal / West Virginia",
  "West_Virginia_Contra" = "Bhopal / West Virginia",
  "Chernobyl"            = "Chernóbil / Fukushima",
  "Fukushima_Contra"     = "Chernóbil / Fukushima",
  "Goiania"              = "Goiânia / Mayapuri",
  "Mayapuri_Contra"      = "Goiânia / Mayapuri",
  "Salisbury"            = "Salisbury"
)

# Nombres legibles para los ejes (los datos usan códigos internos)
etiquetas_casos <- c(
  "Sverdlovsk"           = "Sverdlovsk",
  "Anthrax_Contra"       = "Ántrax EE.UU.",
  "Bhopal"               = "Bhopal",
  "West_Virginia_Contra" = "West Virginia",
  "Chernobyl"            = "Chernóbil",
  "Fukushima_Contra"     = "Fukushima",
  "Goiania"              = "Goiânia",
  "Mayapuri_Contra"      = "Mayapuri",
  "Salisbury"            = "Salisbury"
)

# Paleta por par (una tonalidad por pareja caso/contraparte)
paleta_pares <- c(
  "Sverdlovsk / Ántrax EE.UU." = "#1A5276",  # azul
  "Bhopal / West Virginia"     = "#117A65",  # verde azulado
  "Chernóbil / Fukushima"      = "#B9770E",  # ámbar
  "Goiânia / Mayapuri"         = "#7D3C98",  # púrpura
  "Salisbury"                  = "#922B21"   # granate
)

# Preparar datos para el gráfico
datos_grafico <- datos_todos %>%
  group_by(caso, periodo) %>%
  summarise(total_pub = sum(n_pub, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    periodo = factor(periodo, levels = c("Pre", "Post")),
    caso    = factor(caso, levels = orden_casos),
    par     = factor(grupo_par[as.character(caso)], levels = niveles_par)
  )

# Crear gráfico: color por par, intensidad (alpha) por periodo Pre/Post
grafico <- ggplot(datos_grafico, aes(x = caso, y = total_pub,
                                     fill = par, alpha = periodo)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_text(
    aes(label = total_pub),
    position = position_dodge(width = 0.8),
    vjust    = -0.5,
    size     = 3.2,
    fontface = "bold",
    color    = "grey20",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = paleta_pares, guide = "none") +
  scale_alpha_manual(
    values = c("Pre" = 0.45, "Post" = 1),
    labels = c("Pre (8 años antes)", "Post (8 años después)")
  ) +
  scale_x_discrete(labels = etiquetas_casos) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Análisis Bibliométrico: Publicaciones Pre vs. Post Incidente NRBQ",
    subtitle = "Cada caso principal junto a su contraparte, ordenados por año del incidente",
    x        = "Incidente / Control",
    y        = "Número total de publicaciones",
    alpha    = "Periodo",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5, color = "grey10"),
    plot.subtitle      = element_text(size = 9, hjust = 0.5, color = "grey40"),
    plot.caption       = element_text(size = 8, color = "grey50"),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 35, hjust = 1, size = 9, color = "grey20"),
    axis.text.y        = element_text(size = 9, color = "grey20"),
    axis.title.x       = element_text(size = 10, color = "grey20"),
    axis.title.y       = element_text(size = 10, color = "grey20"),
    legend.position    = "top",
    legend.title       = element_text(size = 9, color = "grey20"),
    legend.text        = element_text(size = 9, color = "grey20")
  )

# Guardar gráfico
archivo_png <- "grafico_bibliometria_TFG.png"
ggsave(archivo_png, plot = grafico, width = 14, height = 7, dpi = 300, bg = "white")
message(paste0("✓ Gráfico guardado en: ", archivo_png))
# ─────────────────────────────────────────────────────────────────────────────
# 7. LOLLIPOP CHART – IRR del Modelo Quasipoisson
# ─────────────────────────────────────────────────────────────────────────────

# Preparar dataframe: mismo orden por pares que el gráfico de barras.
# El eje Y se invierte (rev) para que el primer par quede arriba.
df_lollipop <- resultados |>
  dplyr::filter(!is.na(QuasiP_IRR)) |>
  dplyr::mutate(
    Caso = factor(Caso, levels = rev(orden_casos)),
    par  = factor(grupo_par[as.character(Caso)], levels = niveles_par),
    etiqueta = paste0(
      format(round(QuasiP_IRR, 2), nsmall = 2),
      ifelse(!is.na(QuasiP_pvalor) & QuasiP_pvalor < 0.05, " ***", "")
    )
  )

grafico_irr <- ggplot(df_lollipop, aes(x = QuasiP_IRR, y = Caso)) +
  # Línea de referencia "sin efecto"
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40",
             linewidth = 0.7) +
  # Palo del piruleto (color por par)
  geom_segment(
    aes(x = 1, xend = QuasiP_IRR, y = Caso, yend = Caso, color = par),
    linewidth = 0.8
  ) +
  # Punto: color por par, forma según significación
  geom_point(
    aes(color = par, shape = QuasiP_Significativo),
    size = 4, stroke = 1.1, fill = "white"
  ) +
  # Etiqueta numérica
  geom_text(
    aes(label = etiqueta, color = par),
    hjust  = -0.2,
    size   = 3.2,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_x_log10(
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  scale_y_discrete(labels = etiquetas_casos) +
  scale_color_manual(values = paleta_pares, guide = "none") +
  scale_shape_manual(
    values = c("SÍ" = 16, "NO" = 21),
    labels = c("SÍ" = "Significativo (p < 0.05)", "NO" = "No significativo"),
    name   = "Significación"
  ) +
  labs(
    title    = "Incidence Rate Ratio (IRR) – Modelo Quasipoisson",
    subtitle = "Escala logarítmica | Línea punteada = IRR 1 (sin efecto) | *** p < 0.05",
    x        = "IRR (escala log)",
    y        = NULL,
  ) +
  guides(shape = guide_legend(override.aes = list(color = "grey20"))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5, color = "grey10"),
    plot.subtitle      = element_text(size = 9, hjust = 0.5, color = "grey40"),
    plot.caption       = element_text(size = 8, color = "grey50"),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(size = 9, color = "grey20"),
    axis.title.x       = element_text(size = 10, color = "grey20"),
    legend.position    = "top",
    legend.title       = element_text(size = 9, color = "grey20"),
    legend.text        = element_text(size = 9, color = "grey20")
  )

# Guardar
ggsave("grafico_IRR_Quasipoisson.png", plot = grafico_irr,
       width = 12, height = 5, dpi = 300, bg = "white")
message("✓ Gráfico IRR guardado en: grafico_IRR_Quasipoisson.png")
# ─────────────────────────────────────────────────────────────────────────────
# 8. RESUMEN FINAL EN CONSOLA
# ─────────────────────────────────────────────────────────────────────────────
message("\n")
message("╔══════════════════════════════════════════════════════════╗")
message("║           RESUMEN FINAL DE RESULTADOS                   ║")
message("╚══════════════════════════════════════════════════════════╝")

for (i in seq_len(nrow(resultados))) {
  r <- resultados[i, ]
  message(paste0(
    "\n  ", r$Caso,
    " (", r$Anio_Incidente, ")",
    "  |  Pre: ", r$Pub_Pre, "  Post: ", r$Pub_Post,
    "  |  Binom IRR=", ifelse(is.finite(r$Binom_IRR), round(r$Binom_IRR, 2), "Inf"),
    " p=", format(r$Binom_pvalor, digits = 3, scientific = TRUE),
    " [", r$Binom_Significativo, "]",
    "  |  QP IRR=", ifelse(!is.na(r$QuasiP_IRR), round(r$QuasiP_IRR, 2), "N/A"),
    " p=", ifelse(!is.na(r$QuasiP_pvalor), format(r$QuasiP_pvalor, digits = 3, scientific = TRUE), "N/A"),
    " [", r$QuasiP_Significativo, "]"
  ))
}

message("\n\n════════════════════════════════════════════════════════")
message("  Análisis completado exitosamente.")
message(paste0("  Archivos generados:"))
message(paste0("    • ", archivo_csv))
message(paste0("    • datos_brutos_por_anio_TFG.csv"))
message(paste0("    • ", archivo_png))
message(paste0("    • grafico_IRR_Quasipoisson.png"))
message("════════════════════════════════════════════════════════\n")