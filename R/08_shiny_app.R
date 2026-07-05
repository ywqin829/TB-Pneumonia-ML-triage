# ============================================================
# 08_shiny_app.R
# Interactive Shiny dashboard for TB vs. Pneumonia triage
# Deployed at: https://0pwfel-0-0.shinyapps.io/TB_AI_App/
# Requires: tb_pneumonia_model_V2.rds in app directory
# ============================================================

# ---- Load libraries ----
library(shiny)
library(bslib)
library(xgboost)
library(shapviz)
library(ggplot2)
library(caret)
library(echarts4r)

# ---- Load pre-trained model and dual thresholds ----
model_data <- readRDS("tb_pneumonia_model_V2.rds")
xgb_final_clinical <- model_data$model
lower_thresh <- model_data$lower
upper_thresh <- model_data$upper
xgb_raw_only <- xgb_final_clinical$finalModel

# ---- UI ----
ui <- page_sidebar(
  title = "AI-Powered TB vs. Pneumonia Diagnostic Assistant",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2C3E50"),

  sidebar = sidebar(
    title = "Clinical Input",
    width = 350,
    helpText("Input the 5 parameters as they appear on the routine blood report:"),

    numericInput("RBC", "Red Blood Cell Count (RBC, x10^12/L):", value = 4.5, step = 0.1),
    numericInput("PDW", "Platelet Distribution Width (PDW):", value = 15.0, step = 0.1),
    numericInput("MONO", "Monocyte Percentage (MONO, %):", value = 8.0, step = 0.1),
    numericInput("MPV", "Mean Platelet Volume (MPV, fL):", value = 10.0, step = 0.1),
    numericInput("LYMPH_ABS", "Lymphocyte Count (Absolute, x10^9/L):", value = 1.4, step = 0.1),
    numericInput("MONO_ABS", "Monocyte Count (Absolute, x10^9/L):", value = 0.4, step = 0.05),

    hr(),
    actionButton("predict_btn", "Run AI Diagnosis", class = "btn-danger btn-lg", width = "100%"),

    hr(),
    helpText(HTML("<small style='color:gray;'><b>Disclaimer:</b> For research and
                  educational purposes only. Not for clinical use.</small>"))
  ),

  layout_columns(
    col_widths = c(12, 12),

    # Gauge + risk level card
    card(
      card_header("Diagnostic Prediction", class = "bg-primary text-white"),
      card_body(
        fluidRow(
          column(7, echarts4rOutput("gauge_chart", height = "250px")),
          column(5,
            div(style = "display: flex; flex-direction: column;
                         justify-content: center; height: 250px; padding-left: 10px;",
              h4(textOutput("risk_level"), style = "font-weight: bold;
                 margin-bottom: 20px; line-height: 1.4;"),
              HTML("<div style='font-size: 0.95em; color: #7F8C8D;
                    border-left: 4px solid #BDC3C7; padding-left: 12px;'>
                    <b style='color: #2C3E50;'>Gauge Zones:</b><br/>
                    <span style='color:#E64B35; font-weight:bold;'>● Red</span> :
                      High Risk (TB Suspected)<br/>
                    <span style='color:#95A5A6; font-weight:bold;'>● Gray</span> :
                      Indeterminate (Need Review)<br/>
                    <span style='color:#2980B9; font-weight:bold;'>● Blue</span> :
                      Low Risk (Pneumonia Likely)
                    </div>")
            )
          )
        )
      )
    ),

    # SHAP waterfall card
    card(
      card_header("Explainable AI (XAI) - Individual SHAP Waterfall Plot",
                  class = "bg-dark text-white"),
      card_body(
        h4(textOutput("net_contrib"), style = "text-align: center;
           color: #2C3E50; font-weight: bold; margin-bottom: 5px;"),
        plotOutput("shap_waterfall", height = "450px"),
        helpText(HTML("<div style='text-align: center; font-size: 0.95em; margin-top: 5px;'>
                        <span style='color:#E64B35'><b>Red bars</b></span> = support TB
                        &emsp;|&emsp;
                        <span style='color:#2980B9'><b>Blue bars</b></span> = lean toward Pneumonia
                        &emsp;|&emsp;
                        <b style='color:#2C3E50;'>LMR</b> = LYMPH / MONO
                      </div>"))
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  observeEvent(input$predict_btn, {

    # Compute LMR from lymphocyte and monocyte counts
    safe_mono_abs <- ifelse(input$MONO_ABS == 0, 0.001, input$MONO_ABS)
    calculated_lmr <- input$LYMPH_ABS / safe_mono_abs

    new_patient <- data.frame(
      RBC = input$RBC, PDW = input$PDW, "MONO." = input$MONO / 100,
      MPV = input$MPV, LMR = calculated_lmr, check.names = FALSE
    )

    pred_prob <- predict(xgb_final_clinical, newdata = new_patient,
                         type = "prob")[["TB"]]

    # Risk level text
    output$risk_level <- renderText({
      if (pred_prob < lower_thresh) {
        paste0("LOW RISK (Prob < ", round(lower_thresh * 100, 1),
               "%): TB ruled out. Pneumonia likely.")
      } else if (pred_prob <= upper_thresh) {
        paste0("GRAY ZONE (", round(lower_thresh * 100, 1), "% - ",
               round(upper_thresh * 100, 1),
               "%): Indeterminate. Further tests needed.")
      } else {
        paste0("HIGH RISK (Prob > ", round(upper_thresh * 100, 1),
               "%): TB highly suspected.")
      }
    })

    # Three-color gauge chart
    output$gauge_chart <- renderEcharts4r({
      gauge_color_bands <- list(
        list(lower_thresh, "#2980B9"),   # Blue: low risk
        list(upper_thresh, "#95A5A6"),   # Gray: indeterminate
        list(1, "#E64B35")                # Red: high risk
      )
      e_charts() |>
        e_gauge(
          value = round(pred_prob * 100, 1),
          name = "TB Risk",
          min = 0, max = 100,
          radius = "99%",
          axisLine = list(lineStyle = list(width = 15, color = gauge_color_bands)),
          pointer = list(itemStyle = list(color = "auto")),
          axisLabel = list(distance = 25, fontSize = 10),
          detail = list(
            formatter = htmlwidgets::JS("function(v){return v.toFixed(1)+'%';}"),
            color = "auto", fontSize = 30, fontWeight = "bold",
            offsetCenter = list(0, "65%")
          ),
          title = list(offsetCenter = list(0, "30%"),
                       fontWeight = "bold", fontSize = 18, color = "auto")
        )
    })

    # SHAP waterfall
    patient_mat <- as.matrix(new_patient)
    sv_new <- shapviz(xgb_raw_only, X_pred = patient_mat, X = patient_mat)
    sv_new$S <- -sv_new$S
    sv_new$baseline <- -sv_new$baseline

    net_shap <- sum(sv_new$S)
    output$net_contrib <- renderText({
      if (net_shap > 0) {
        paste0("Net SHAP = +", round(net_shap, 3),
               " — Strongly points to TB")
      } else {
        paste0("Net SHAP = ", round(net_shap, 3),
               " — Points to Pneumonia")
      }
    })

    output$shap_waterfall <- renderPlot({
      suppressWarnings({
        mono_idx <- which(colnames(sv_new$S) == "MONO.")
        colnames(sv_new$S)[mono_idx] <- "MONO%"
        colnames(sv_new$X)[mono_idx] <- "MONO%"

        sv_waterfall(sv_new, fill_colors = c("#E64B35", "#2980B9")) +
          theme_bw() +
          theme(
            plot.title = element_text(face = "bold", size = 16,
                                      hjust = 0.5, margin = margin(b = 10)),
            axis.text.y = element_text(face = "bold", size = 14, color = "black"),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.title.x = element_text(face = "bold", size = 14)
          ) +
          labs(title = "Physiological Drivers for Current Patient",
               x = "SHAP Value (Impact on TB Probability)")
      })
    })
  })
}

shinyApp(ui, server)
