library(readxl)
library(ggplot2)
library(shiny)
library(shinyWidgets)
library(bslib)
library(shinyjs)
library(lubridate)
library(dplyr)



procesar_datos <- function(data, variables_interes) {
  library(dplyr)
  
  # Eliminar valores específicos no válidos
  if ("EDAD_DIAGNOSTICO" %in% colnames(data)) {
    data <- data %>% filter(EDAD_DIAGNOSTICO != 999)
  }
  if ("Nx" %in% colnames(data)) {
    data <- data %>% filter(Nx != 999)
  }
  if ("TAMAÑO" %in% colnames(data)) {
    data <- data %>% filter(TAMAÑO != 999)
  }
  
  # Imputación de valores perdidos y "Desconocido"
  for (var in variables_interes) {
    if (var %in% colnames(data)) {
      if (is.numeric(data[[var]])) {
        # Imputar NA numéricos con la media
        data[[var]][is.na(data[[var]])] <- mean(data[[var]], na.rm = TRUE)
      } else if (is.factor(data[[var]]) || is.character(data[[var]])) {
        # Reemplazar "Desconocido" con NA
        data[[var]][data[[var]] == "Desconocido"] <- NA
        
        # Calcular la moda e imputar
        moda <- names(which.max(table(data[[var]])))
        data[[var]][is.na(data[[var]])] <- moda
      }
    }
  }
  
  # Detección y eliminación de valores atípicos (solo en numéricos)
  for (var in variables_interes) {
    if (var %in% colnames(data) && is.numeric(data[[var]])) {
      Q1 <- quantile(data[[var]], 0.25, na.rm = TRUE)
      Q3 <- quantile(data[[var]], 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      limite_inferior <- Q1 - 1.5 * IQR
      limite_superior <- Q3 + 1.5 * IQR
      data <- data %>% filter(data[[var]] >= limite_inferior & data[[var]] <= limite_superior)
    }
  }
  
  # Seleccionar solo las variables de interés
  data <- dplyr::select(data, dplyr::all_of(variables_interes))
  
  # Eliminar niveles no utilizados en factores
  data <- droplevels(data)
  
  return(data)
}

# Función para preprocesar el dataset
preprocesar_datos <- function(data) {
  
  # Verificar si las columnas existen antes de modificarlas
  if ("GRADO" %in% colnames(data)) {
    data <- data %>%
      mutate(GRADO_CAT = case_when(
        GRADO %in% c(1, 2) ~ "1/2",
        GRADO == 3 ~ "3",
        TRUE ~ NA_character_
      )) %>%
      mutate(GRADO_CAT = as.factor(GRADO_CAT))
  }
  
  if ("HISTOLOGIA" %in% colnames(data)) {
    data <- data %>%
      mutate(HISTOLOGIA_CAT = case_when(
        HISTOLOGIA == "DUCTAL" ~ "DUCTAL",
        TRUE ~ "NO DUCTAL"
      )) %>%
      mutate(HISTOLOGIA_CAT = as.factor(HISTOLOGIA_CAT))
  }
  
  if ("Ki67" %in% colnames(data)) {
    data <- data %>%
      mutate(Ki67_CAT = case_when(
        Ki67 <= 14 ~ "<=14",
        Ki67 > 14 ~ ">14",
        TRUE ~ NA_character_
      )) %>%
      mutate(Ki67_CAT = as.factor(Ki67_CAT))
  }
  
  if ("FENOTIPO_IHQ" %in% colnames(data)) {
    data <- data %>%
      mutate(FENOTIPO_IHQ_CAT = case_when(
        FENOTIPO_IHQ %in% c("LUMINAL A", "LUMINAL B") ~ "LUMINAL",
        FENOTIPO_IHQ %in% c("LUMINAL-HER2", "HER2") ~ "HER2",
        TRUE ~ "TRIPLE NEGATIVO"
      )) %>%
      mutate(FENOTIPO_IHQ_CAT = as.factor(FENOTIPO_IHQ_CAT))
  }
  
  if ("ESTADIO" %in% colnames(data)) {
    data <- data %>%
      mutate(ESTADIO_CAT = case_when(
        ESTADIO %in% c("I", "IA", "IB") ~ "I",
        ESTADIO %in% c("II", "IIA", "IIB") ~ "II",
        ESTADIO %in% c("III", "IIIA", "IIIB", "IIIC") ~ "III",
        TRUE ~ NA_character_
      )) %>%
      mutate(ESTADIO_CAT = as.factor(ESTADIO_CAT))
  }
  
  if ("TIPO_QT_NEOADY_OK" %in% colnames(data)) {
    data <- data %>%
      mutate(QT_CAT = case_when(
        TIPO_QT_NEOADY_OK %in% c("ANTRACICLINAS", "ANTRACICLINAS-TAXANOS", "ANTRACICLINAS-TAXANOS-PLATINOS", "TAXANOS") ~ "NO antiHER2",
        TIPO_QT_NEOADY_OK %in% c("QUIMIOTERAPIA-ANTIHER2", "QUIMIOTERAPIA-ANTIHER2DOBLE") ~ "antiHER2",
        TIPO_QT_NEOADY_OK %in% c("QUIMIOTERAPIA-ANTIHER2-INMUNOTERAPIA", "QUIMIOTERAPIA-INMUNOTERAPIA") ~ "INMUNO",
        TRUE ~ NA_character_
      )) %>%
      mutate(QT_CAT = as.factor(QT_CAT))
  }
  
  if ("No_ANTIBIOTICOS" %in% colnames(data)) {
    data <- data %>%
      mutate(No_ANTIBIOTICOS_CAT = case_when(
        No_ANTIBIOTICOS == 1 ~ "Una vez",
        No_ANTIBIOTICOS >= 2 ~ "Dos veces o más",
        TRUE ~ NA_character_
      )) %>%
      mutate(No_ANTIBIOTICOS_CAT = as.factor(No_ANTIBIOTICOS_CAT)) %>%
      mutate(No_ANTIBIOTICOS_CAT = droplevels(No_ANTIBIOTICOS_CAT))
  }
  
  if ("RCB" %in% colnames(data)) {
    data <- data %>%
      mutate(RCP = case_when(
        RCB %in% c("RCB-0", "RCB-I") ~ "SI",
        RCB %in% c("RCB-II", "RCB-III") ~ "NO",
        TRUE ~ NA_character_
      )) %>%
      mutate(RCP = as.factor(RCP))
  }
  
  return(data)
}


muestra_700_UGC_perdidos <- read_excel("muestra_700_UGC_perdidos.xlsx")

# Convertir variables character en factor
muestra_700_UGC_perdidos <- muestra_700_UGC_perdidos %>%
  mutate_if(is.character, as.factor)

muestra_700_UGC_perdidos$Tn <- factor(muestra_700_UGC_perdidos$Tn)
muestra_700_UGC_perdidos$GRADO <- factor(muestra_700_UGC_perdidos$GRADO)
muestra_700_UGC_perdidos$Nx <- factor(muestra_700_UGC_perdidos$Nx)

variables_a_analizar <- c("EDAD_DIAGNOSTICO", "ESTADO_MENOPAUSICO", "HISTOLOGIA", "GRADO", "KI67_DICOT", "Ki67", "FENOTIPO_IHQ", "TAMAÑO", "Tn", "Nx", "ESTADIO", "TIPO_CIRUGIA", "TIPO_QT_NEOADY_OK", "RCB", "HT_ADY_OK", "RE", "RP", "RETRASOS", "No_ANTIBIOTICOS", "FECHA_ULTIMO_CONTROL", "FECHA_INICIO_QT", "FECHA_DIAGNOSTICO", "ESTADO_ULTIMO_\nCONTROL", "ANTIBIOTICO")

df_procesado <- procesar_datos(muestra_700_UGC_perdidos, variables_a_analizar)

df_procesado$EDAD_DIAGNOSTICO <- cut(df_procesado$EDAD_DIAGNOSTICO, breaks = c(-Inf, 39, 49, 59,69,79, Inf), labels = c("0-39", "40-49", "50-59", "60-69","70-79","≥80"))

df_procesado$EDAD_DIAGNOSTICO <- as.factor(df_procesado$EDAD_DIAGNOSTICO)

df_procesado <- df_procesado %>%
  mutate(
    FECHA_DIAGNOSTICO = dmy(FECHA_DIAGNOSTICO),
    FECHA_DIAGNOSTICO = cut(
      year(FECHA_DIAGNOSTICO),
      breaks = seq(2000, 2030, by = 10),
      labels = c("2000-2009", "2010-2019", "2020-2029"),
      right = FALSE
    ),
    FECHA_DIAGNOSTICO = as.factor(FECHA_DIAGNOSTICO)
  )

df_procesado <- df_procesado %>%
  mutate(
    FECHA_ULTIMO_CONTROL = dmy(FECHA_ULTIMO_CONTROL),
    FECHA_ULTIMO_CONTROL = cut(
      year(FECHA_ULTIMO_CONTROL),
      breaks = seq(2000, 2030, by = 10),
      labels = c("2000-2009", "2010-2019", "2020-2029"),
      right = FALSE
    ),
    FECHA_ULTIMO_CONTROL = as.factor(FECHA_ULTIMO_CONTROL)
  )

df_procesado <- df_procesado %>%
  mutate(
    FECHA_INICIO_QT = dmy(FECHA_INICIO_QT),
    FECHA_INICIO_QT = cut(
      year(FECHA_INICIO_QT),
      breaks = seq(2000, 2030, by = 10),
      labels = c("2000-2009", "2010-2019", "2020-2029"),
      right = FALSE
    ),
    FECHA_INICIO_QT = as.factor(FECHA_INICIO_QT)
  )

# Seleccionar variables numéricas
numeric_vars <- sapply(df_procesado, is.numeric)

# Excluir 'Ki67' y 'No_ANTIBIOTICOS' de la normalización
vars_a_normalizar <- setdiff(names(df_procesado)[numeric_vars], c("Ki67", "No_ANTIBIOTICOS"))

dataset_normalizar_no <- df_procesado

media_tam <- mean(dataset_normalizar_no$TAMAÑO)
sd_tam <- sd(dataset_normalizar_no$TAMAÑO)

media_re <- mean(dataset_normalizar_no$RE)
sd_re <- sd(dataset_normalizar_no$RE)

# Normalizar solo esas variables
normalized_data <- scale(df_procesado[, vars_a_normalizar])

# Reemplazar las columnas originales por las normalizadas
df_procesado[, vars_a_normalizar] <- normalized_data

dataset_recodificado <- preprocesar_datos(df_procesado)

# Eliminar niveles no usados
dataset_recodificado <- droplevels(dataset_recodificado)

# Eliminar columnas innecesarias
dataset_recodificado <- subset(dataset_recodificado, select = -c(ANTIBIOTICO, GRADO, HISTOLOGIA, Ki67, KI67_DICOT, FENOTIPO_IHQ, ESTADIO, TIPO_QT_NEOADY_OK,No_ANTIBIOTICOS, RCB))

# Arreglar nombre de columna malformado
names(dataset_recodificado)[names(dataset_recodificado) == "ESTADO_ULTIMO_\nCONTROL"] <- "ESTADO_ULTIMO_CONTROL"

dataset_modelo <- na.omit(dataset_recodificado)
dataset_modelo <- droplevels(dataset_modelo)

library(shiny)

# Entrenar modelo con fórmula dada
modelo_rlog <- glm(RCP ~ EDAD_DIAGNOSTICO + TAMAÑO + Tn + RE + ESTADO_ULTIMO_CONTROL + QT_CAT,
                   data = dataset_modelo, family = binomial)

# Temas
tema_claro <- bs_theme(version = 5, bootswatch = "flatly", primary = "#007BFF")
tema_oscuro <- bs_theme(version = 5, bootswatch = "darkly", primary = "#00C2FF")

ui <- fluidPage(
  useShinyjs(),
  theme = tema_claro,
  
  tags$head(
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css"),
    tags$link(rel = "stylesheet", href = "https://www.w3schools.com/w3css/4/w3.css"),
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Open+Sans&display=swap"),
    
    tags$style(HTML("
      body {
        font-family: 'Open Sans', sans-serif;
      }
      body, .well, .panel, .form-control, .btn, .custom-card {
        transition: background-color 0.4s ease, color 0.4s ease, border-color 0.4s ease;
      }

      .dark-theme {
        background-color: #121212 !important;
        color: #ffffff !important;
      }

      .dark-theme a, 
      .dark-theme .form-control, 
      .dark-theme .btn {
        color: #ffffff !important;
        background-color: #1e1e1e !important;
        border-color: #444 !important;
      }

      .dark-theme .custom-card {
        background-color: #1e1e1e !important;
      }
      
      .grafico-redondeado img {
        border-radius: 10px;
        border: 1px solid #ccc;
        box-shadow: 0 2px 6px rgba(0, 0, 0, 0.1); /* opcional para resaltar más */
      }

      .theme-toggle-btn {
        position: fixed;
        top: 15px;
        right: 20px;
        z-index: 9999;
        font-size: 20px;
        background-color: transparent;
        color: #007BFF;
        border: none;
      }

      .theme-toggle-btn:hover {
        color: #0056b3;
      }
    ")),
    
    tags$script(HTML("
      function retheme() {
        const body = document.body;
        const isDark = body.classList.toggle('dark-theme');
        Shiny.setInputValue('theme_toggle', isDark ? 'oscuro' : 'claro', {priority: 'event'});
      }

      $(document).on('shown.bs.modal', function() {
        if (document.body.classList.contains('dark-theme')) {
          $('.modal-content').css({
            'background-color': '#1e1e1e',
            'color': '#ffffff',
            'border': '1px solid #555'
          });
          $('.modal-header, .modal-footer').css({
            'background-color': '#1e1e1e',
            'border-color': '#444'
          });
        } else {
          $('.modal-content, .modal-header, .modal-footer').removeAttr('style');
        }
      });
    "))
  ),
  
  
  # Botón de cambio de tema (arriba derecha)
  tags$button(
    class = "theme-toggle-btn fa fa-adjust",
    title = "Cambiar tema",
    onclick = "retheme()"
  ),
  
  br(), br(),
  
  titlePanel(div(
    tags$h2("Predicción de Respuesta Completa Patológica (RCP)", style = "text-align:center; font-weight:700; color:#007BFF")
  )),
  
  br(),
  
  # Ayuda centrada
  div(class = "text-center",
      actionLink("ayuda", "¿Necesitas ayuda para entender las variables?")
  ),
  
  br(),
  
  # Inputs del modelo
  fluidRow(
    column(6,
           pickerInput("edad", "Edad al diagnóstico", 
                       choices = sort(unique(dataset_modelo$EDAD_DIAGNOSTICO)),
                       options = list(style = "btn-info")),
           sliderInput("tamano", "Tamaño del tumor (cm)", 
                       min = min(dataset_normalizar_no$TAMAÑO),
                       max = max(dataset_normalizar_no$TAMAÑO),
                       value = 1, step = 0.1),
           selectInput("tn", "Clasificación del tumor primario", choices = sort(unique(dataset_modelo$Tn)))
    ),
    column(6,
           div(style = "display: flex; flex-direction: column; align-items: flex-end;",
               sliderInput("re", "Expresión de RE", 
                           min = min(dataset_normalizar_no$RE),
                           max = max(dataset_normalizar_no$RE),
                           value = 0, step = 1),
               selectInput("estado", "Último estado clínico", choices = sort(unique(dataset_modelo$ESTADO_ULTIMO_CONTROL))),
               selectInput("qt", "Tipo de Quimioterapia", choices = sort(unique(dataset_modelo$QT_CAT)))
           )
    )
  )
  ,
  
  br(),
  
  # Botón de predicción
  fluidRow(
    column(12, align = "center",
           actionButton("predecir", label = "Predecir", class = "btn btn-success btn-lg")
    )
  ),
  
  br(),
  
  # Resultados
  fluidRow(
    column(6,
           uiOutput("card_resultado")
    ),
    column(6,div(class = "grafico-redondeado",
                 plotOutput("grafico_prob", height = "250px")
    )
    )
  )
)


server <- function(input, output, session) {
  
  observeEvent(input$theme_toggle, {
    # Puedes hacer algo con el valor, por ejemplo:
    cat("Tema actual:", input$theme_toggle, "\n")
  })
  
  # Mostrarlo en la UI
  output$theme_status <- renderText({
    paste("Tema actual:", input$theme_toggle)
  })
  
  observeEvent(input$ayuda, {
    showModal(modalDialog(
      title = "Ayuda sobre las variables",
      tags$ul(
        tags$li("Edad al diagnóstico: Rango de edad al momento del diagnóstico."),
        tags$br(),
        tags$li("Tamaño: Tamaño del tumor en cm."),
        tags$br(),
        tags$li("Tn: Clasificación del tumor primario."),
        tags$br(),
        tags$li("RE: Nivel de expresión de receptores de estrógeno."),
        tags$br(),
        tags$li("Estado último control: Último estado clínico registrado."),
        tags$br(),
        tags$li("QT_CAT: Tipo de quimioterapia neoadyuvante recibida.")
      ),
      easyClose = TRUE,
      footer = modalButton("Cerrar")
    ))
  })
  
  prediccion <- eventReactive(input$predecir, {
    tamaño_input <- (input$tamano - media_tam) / sd_tam
    re_input <- (input$re - media_re) / sd_re
    
    nuevo <- data.frame(
      EDAD_DIAGNOSTICO = factor(input$edad, levels = levels(dataset_modelo$EDAD_DIAGNOSTICO)),
      TAMAÑO = tamaño_input,
      Tn = factor(input$tn, levels = levels(dataset_modelo$Tn)),
      RE = re_input,
      ESTADO_ULTIMO_CONTROL = factor(input$estado, levels = levels(dataset_modelo$ESTADO_ULTIMO_CONTROL)),
      QT_CAT = factor(input$qt, levels = levels(dataset_modelo$QT_CAT))
    )
    
    prob <- predict(modelo_rlog, newdata = nuevo, type = "response")
    clase <- ifelse(prob >= 0.5, "Sí", "No")
    
    list(probabilidad = prob, clase = clase)
  })
  
  output$card_resultado <- renderUI({
    req(prediccion())
    
    prob <- round(prediccion()$probabilidad * 100, 1)
    clase <- prediccion()$clase
    
    color <- if (clase == "Sí") "#00AEDA" else "#dc3545"
    mensaje <- if (clase == "Sí") "Probabilidad alta de RCP" else "Probabilidad baja de RCP"
    confianza <- ifelse(prob > 80, "Alta confianza", ifelse(prob > 60, "Moderada", "Baja"))
    
    div(
      style = paste0("background-color:", color, "; color: white; padding: 25px; border-radius: 10px; text-align: center;"),
      tags$strong(h3(mensaje)),
      tags$p(paste("Probabilidad estimada:", prob, "%")),
      tags$strong(paste("Nivel de confianza del modelo:", confianza)),
      tags$hr(),
      tags$small("Esta estimación se basa en un modelo clínico de regresión logística.")
    )
  })
  
  output$grafico_prob <- renderPlot({
    req(prediccion())
    prob <- prediccion()$probabilidad
    
    barplot(
      height = c(prob, 1 - prob),
      names.arg = c("RCP: Sí", "RCP: No"),
      col = c("#00AEDA", "#002859"),
      ylim = c(0, 1),
      main = "Distribución de la Probabilidad"
    )
  })
}

shinyApp(ui = ui, server = server)
