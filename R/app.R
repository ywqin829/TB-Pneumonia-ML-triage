## Local Shiny research calculator for the signed final three-feature model.
## Start with shiny::runApp("D:/科研/6投稿/分诊模型/rerun", launch.browser=TRUE).
suppressPackageStartupMessages(library(shiny))
app_root <- getwd()
source(file.path(app_root, "app_predict.R"), local = TRUE)
bundle <- load_tb_bundle(app_root)

ui <- fluidPage(
  tags$head(tags$style(HTML(paste(
    "body {background:#f6f7f9;color:#172334;font-family:Arial,'Microsoft YaHei',sans-serif;}",
    ".container-fluid {max-width:1050px;padding-top:22px;} .well {background:white;}",
    ".result {padding:18px;border-radius:8px;background:white;border-left:7px solid;}",
    ".score {font-size:36px;font-weight:700;} .note {color:#516174;font-size:14px;}",
    ".caution {padding:12px;background:#fff4dc;border-radius:6px;}",
    "h2 {font-size:25px;} h4 {font-weight:600;}")))),
  titlePanel("结核/肺炎重跑模型 · TB/Pneumonia Research Calculator"),
  p(class = "note", "2026-10-04 交付 / Delivery · 3-feature XGBoost · 平衡分类方案 / Balanced classification"),
  sidebarLayout(
    sidebarPanel(
      numericInput("mcv", "MCV（与源数据相同量纲 / source scale）", value = 90, min = 1, step = 0.1),
      numericInput("pdw", "PDW（与源数据相同量纲 / source scale）", value = 12, min = 0.01, step = 0.1),
      numericInput("mono_percent", "单核细胞百分比 / MONO (%)", value = 6, min = 0, max = 100, step = 0.1),
      actionButton("calculate", "计算 / Calculate", class = "btn-primary"),
      hr(),
      p(class = "note", "MONO% 输入 6 表示 6%，模型接收 0.06。MCV/PDW 必须与建模数据量纲一致；CSV 未提供仪器和单位元数据。"),
      p(class = "note", "MONO 6 means 6%; it is converted to 0.06. Match MCV/PDW to the source scale. Analyzer and unit metadata were not supplied."),
      h4("训练数据范围 / Training ranges"),
      tableOutput("ranges"), width = 4
    ),
    mainPanel(
      uiOutput("prediction"),
      br(),
      div(class = "caution",
          strong("结果使用范围 / Interpretation"),
          p("这是内部重分析模型输出，不是确诊结论，也不是已校准的个体患病风险。测试集显示预测概率整体偏高；低分区仍有结核病例。"),
          p("This internally reanalyzed score is not a diagnosis or a calibrated individual risk. Test-set calibration showed probability overestimation, and TB cases remained in the low-score zone.")),
      h4("当前平衡二分类规则 / Current balanced classification"),
      p(sprintf("验证集Youden阈值 / Validation Youden threshold: %.6f", bundle$thresholds[["balanced"]])),
      p("达到阈值时模型倾向结核，低于阈值时模型倾向肺炎。测试敏感度84.6%、特异度71.0%；这不是确诊结论。"),
      p(class = "note", "At or above the threshold the model favors TB; below it the model favors pneumonia. Test sensitivity 84.6%, specificity 71.0%; this is not a diagnosis."),
      p(class = "note", "本次目标调整发生在已查看测试结果之后，虽然Youden阈值取自验证集，也不能宣称为事前预定的主方案。 / The primary objective was changed after test results were reviewed; the validation-derived Youden threshold does not make this a prospectively specified primary scheme."),
      h4("补充分诊分区 / Supplementary triage zones"),
      tableOutput("zone_rules"),
      p(class = "note", "三区边界沿用验证集敏感度≥90%与特异度≥90%阈值，仅作补充。低阈值测试敏感度91.6%、特异度50.7%；高阈值测试特异度85.7%。"),
      p(class = "note", "Validation targets do not guarantee the same sensitivity or specificity on a new dataset. Further diagnostic assessment is required in all zones."),
      p(class = "note", "同一来源队列已用于原研究。当前测试划分属于内部重分析，不能称为独立外部或时间验证。 / The source cohort was used in the original study; this split is internal reanalysis, not independent external or temporal validation."),
      width = 8
    )
  )
)

server <- function(input, output, session) {
  result <- eventReactive(input$calculate, {
    validate(need(!is.null(input$mcv) && !is.null(input$pdw) && !is.null(input$mono_percent),
                  "请填写全部数值 / Complete all inputs."))
    tryCatch(predict_tb_inputs(bundle, input$mcv, input$pdw, input$mono_percent),
             error = function(e) list(error = conditionMessage(e)))
  }, ignoreInit = FALSE)
  output$prediction <- renderUI({
    if (input$calculate == 0) return(p("输入血常规指标后点击计算。 / Enter the CBC values and calculate."))
    prediction <- result()
    if (!is.null(prediction$error)) return(div(class = "caution", prediction$error))
    styles <- c(TB = "#e34a33", Pneumonia = "#377eb8")
    classes <- c(TB = "模型倾向结核 / Model favors TB", Pneumonia = "模型倾向肺炎 / Model favors pneumonia")
    labels <- c(low = "低分区 / Low-score zone", intermediate = "中间区 / Intermediate zone",
                 high = "高分区 / High-score zone")
    div(class = "result", style = paste0("border-left-color:", styles[[prediction$model_class]], ";"),
        h4(classes[[prediction$model_class]]),
        div(class = "score", sprintf("%.1f%%", prediction$score * 100)),
        p("模型输出 / Model output score"),
        p(sprintf("平衡分类阈值 / Balanced threshold: %.6f", bundle$thresholds[["balanced"]])),
        p(paste("补充分诊分区 / Supplementary triage:", labels[[prediction$zone]])),
        if (length(prediction$outside_training_range)) {
          p(class = "caution", paste("超出训练数据范围 / Outside training range:",
                                      paste(prediction$outside_training_range, collapse = ", ")))
        },
        p(class = "note", "请结合临床、影像和病原学检查。 / Interpret with clinical, radiological, and microbiological assessment."))
  })
  output$ranges <- renderTable({
    data.frame(Input = c("MCV", "PDW", "MONO (%)"),
               Min = round(bundle$ranges[, 1], 3), Max = round(bundle$ranges[, 2], 3),
               row.names = NULL)
  })
  output$zone_rules <- renderTable({
    data.frame(Zone = c("低 / Low", "中间 / Intermediate", "高 / High"),
               Score = c(sprintf("p < %.6f", bundle$thresholds[["low"]]),
                 sprintf("%.6f ≤ p < %.6f", bundle$thresholds[["low"]], bundle$thresholds[["high"]]),
                 sprintf("p ≥ %.6f", bundle$thresholds[["high"]])),
               Test_TB = c("12 / 182", "45 / 162", "86 / 134"))
  }, striped = TRUE)
}

shinyApp(ui, server)
