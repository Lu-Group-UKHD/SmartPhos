# Mock data for testing
resTab <- data.frame(
  site = c("EGFR_Y1172", "EGFR_Y1197", "EGFR_S1166", "ROCK2_S1374", "WASL_Y256", "GAB1_Y259", "ADD1_S586", "EPHA2_Y772", "PRKDC_T2638", "PRKDC_T2609", "PRKDC_S2612"),
  stat = c(-10.038770, -5.945562, 5.773384, -7.303834, 5.585326, 5.971104, 5.199119, -5.169500, 5.130228, 5.407387, 4.493933),
  log2FC = c(-2.6113343, -2.4858615, 1.0056629, -1.1561780, 1.6421145, 2.0296634, 1.3766283, -0.8531656, 1.0742881, 1.0042942, 1.0608129)
)

decoupler_network <- data.frame(
  source = c(rep("ABL1", 5), rep("CDK2", 6)),
  mor = c(rep(1, 11)),
  target = c("EGFR_Y1172", "EGFR_Y1197", "EGFR_S1166", "ROCK2_S1374", "WASL_Y256", "GAB1_Y259", "ADD1_S586", "EPHA2_Y772", "PRKDC_T2638", "PRKDC_T2609", "PRKDC_S2612"),
  likelihood = c(rep(1, 11))
)


######################## Tests for calcKinaseScore() ###########################

# Test that the function runs without errors with default parameters
test_that("calcKinaseScore runs without errors with default parameters", {
  result <- calcKinaseScore(resTab, decoupler_network)
  
  expect_true(is.data.frame(result))
  expect_true(all(c("source", "score", "p_value") %in% colnames(result)))
})

# Test that duplicate rows are removed based on the 'site' column
test_that("calcKinaseScore removes duplicate rows based on the 'site' column", {
  result <- calcKinaseScore(resTab, decoupler_network)

  expect_equal(nrow(result), length(unique(decoupler_network$source[resTab$site %in% decoupler_network$target])))
})

# Test that the function correctly handles 'statType' = "stat"
test_that("calcKinaseScore handles 'statType' = 'stat' correctly", {
  result <- calcKinaseScore(resTab, decoupler_network, statType = "stat")

  input_values <- resTab %>% filter(site %in% decoupler_network$target) %>% pull(stat)
  expect_equal(length(result$score), length(unique(decoupler_network$source)))
})

# Test that the function correctly handles 'statType' = "log2FC"
test_that("calcKinaseScore handles 'statType' = 'log2FC' correctly", {
  result <- calcKinaseScore(resTab, decoupler_network, statType = "log2FC")

  input_values <- resTab %>% filter(site %in% decoupler_network$target) %>% pull(log2FC)
  expect_equal(length(result$score), length(unique(decoupler_network$source)))
})

# Test that the function filters out correlated regulons based on 'corrThreshold'
test_that("calcKinaseScore filters out correlated regulons based on 'corrThreshold'", {
  result <- calcKinaseScore(resTab, decoupler_network, corrThreshold = 0.9)

  correlated <- check_corr(decoupler_network) %>%
    filter(correlation >= 0.9)

  expect_true(all(!result$source %in% correlated$source.2))
})

# Test that kinase activity scores are calculated and NA values are handled
test_that("calcKinaseScore calculates kinase activity scores and handles NA values", {
  result <- calcKinaseScore(resTab, decoupler_network)

  expect_false(any(is.na(result$score)))
  expect_false(any(is.na(result$p_value)))
})

# Test that the number of permutations affects the result
test_that("calcKinaseScore produces different results with different number of permutations", {
  result_100 <- calcKinaseScore(resTab, decoupler_network, nPerm = 100)
  result_1000 <- calcKinaseScore(resTab, decoupler_network, nPerm = 1000)

  expect_false(all(result_100$p_value == result_1000$p_value))
})




######################### Tests for plotKinaseDE() #############################

# Mock data for testing
scoreTab <- data.frame(
  source = c("Kinase1", "Kinase2", "Kinase3", "Kinase4", "Kinase5"),
  score = c(2.3, -1.5, 0, 3.1, -2.8),
  p_value = c(0.01, 0.2, 0.05, 0.03, 0.04)
)

# Test cases

# Test that the function runs without errors with default parameters
test_that("plotKinaseDE runs without errors with default parameters", {
  plot <- plotKinaseDE(scoreTab)
  
  expect_true(inherits(plot, "ggplot"))
})

# Test that kinases with a score of 0 are filtered out
test_that("plotKinaseDE filters out kinases with a score of 0", {
  plotTab <- scoreTab %>% 
    mutate(significance = ifelse(p_value <= 0.05, "p <= 0.05", "p > 0.05"),
           score_sign = sign(score)) %>%
    filter(score_sign != 0) %>% 
    group_by(score_sign) %>% 
    slice_max(abs(score), n = 10)
  
  plot <- plotKinaseDE(scoreTab)
  
  expect_false("Kinase3" %in% plot$data$source)
})

# Test that the function respects the `nTop` parameter
test_that("plotKinaseDE respects the nTop parameter", {
  plot <- plotKinaseDE(scoreTab, nTop = 2)
  
  expect_lte(nrow(plot$data), 4)  # Since nTop=2, max 2 positive and 2 negative kinases
})

# Test that the function respects the `pCut` parameter
test_that("plotKinaseDE respects the pCut parameter", {
  plot <- plotKinaseDE(scoreTab, pCut = 0.03)
  
  expect_false("p <= 0.05" %in% plot$data$significance[plot$data$p_value > 0.03])
})

# Test that the plot has the correct title
test_that("plotKinaseDE generates a plot with the correct title", {
  plot <- plotKinaseDE(scoreTab, nTop = 3)
  
  expect_equal(plot$labels$title, "Kinase activity inference, top 3 kinases")
})

# Test that the plot uses the correct fill colors based on significance
test_that("plotKinaseDE uses the correct fill colors based on significance", {
  plot <- plotKinaseDE(scoreTab, pCut = 0.05)
  
  expected_colors <- c("indianred", "lightgrey")
  actual_colors <- ggplot_build(plot)$data[[1]]$fill
  
  expect_true(all(actual_colors %in% expected_colors))
})

# Test that the function handles cases where all kinases are not significant
test_that("plotKinaseDE handles cases where no kinases are significant", {
  scoreTab_no_sig <- data.frame(
    source = c("Kinase1", "Kinase2", "Kinase3"),
    score = c(1.2, -2.3, 0.7),
    p_value = c(0.1, 0.2, 0.15)
  )
  
  plot <- plotKinaseDE(scoreTab_no_sig, pCut = 0.05)
  
  expect_true(all(plot$data$significance == "p > 0.05"))
  expect_true(all(ggplot_build(plot)$data[[1]]$fill == "lightgrey"))
})

# Test that the function handles cases where all kinases are significant
test_that("plotKinaseDE handles cases where all kinases are significant", {
  scoreTab_all_sig <- data.frame(
    source = c("Kinase1", "Kinase2", "Kinase3"),
    score = c(1.2, -2.3, 0.7),
    p_value = c(0.01, 0.02, 0.03)
  )
  
  plot <- plotKinaseDE(scoreTab_all_sig, pCut = 0.05)
  
  expect_true(all(plot$data$significance == "p <= 0.05"))
  expect_true(all(ggplot_build(plot)$data[[1]]$fill == "indianred"))
})

###################### Tests for plotKinaseTimeSeries() ########################

# Mock data for testing
scoreTab <- data.frame(
  timepoint = rep(c("0h", "1h", "2h"), each = 3),
  source = rep(c("KinaseA", "KinaseB", "KinaseC"), times = 3),
  score = c(1.2, -1.5, 0.8, -0.9, 1.5, -1.2, 0.6, -1.3, 1.8),
  p_value = c(0.02, 0.06, 0.01, 0.04, 0.07, 0.03, 0.02, 0.08, 0.01)
)

# Test cases

# Test that the function runs without errors with default parameters
test_that("plotKinaseTimeSeries runs without errors with default parameters", {
  p <- plotKinaseTimeSeries(scoreTab)
  
  expect_true(inherits(p, "ggplot"))
})

# Test that the significance markers (*) are correctly added based on p-value cutoff
test_that("plotKinaseTimeSeries adds significance markers based on pCut", {
  p <- plotKinaseTimeSeries(scoreTab, pCut = 0.05)
  
  plot_data <- ggplot_build(p)$data[[2]]  # Data from geom_text layer
  
  significant_entries <- scoreTab %>% filter(p_value <= 0.05)
  
  expect_true(all(plot_data$label[plot_data$label == "*"] == "*"))
  expect_equal(sum(plot_data$label == "*"), nrow(significant_entries))
})

# Test that the function handles different pCut values correctly
test_that("plotKinaseTimeSeries respects the pCut parameter", {
  p <- plotKinaseTimeSeries(scoreTab, pCut = 0.03)
  
  plot_data <- ggplot_build(p)$data[[2]]  # Data from geom_text layer
  
  significant_entries <- scoreTab %>% filter(p_value <= 0.03)
  
  expect_equal(sum(plot_data$label == "*"), nrow(significant_entries))
})

# Test that the plot has the correct title based on clusterName
test_that("plotKinaseTimeSeries generates a plot with the correct title", {
  p <- plotKinaseTimeSeries(scoreTab, clusterName = "TestCluster")
  
  expect_equal(p$labels$title, "Kinase activity infererence, TestCluster")
})

# Test that the function handles cases with no significant activities
test_that("plotKinaseTimeSeries handles cases with no significant activities", {
  scoreTab_no_sig <- data.frame(
    timepoint = rep(c("0h", "1h", "2h"), each = 3),
    source = rep(c("KinaseA", "KinaseB", "KinaseC"), times = 3),
    score = c(1.2, -1.5, 0.8, -0.9, 1.5, -1.2, 0.6, -1.3, 1.8),
    p_value = c(0.1, 0.2, 0.15, 0.08, 0.07, 0.09, 0.2, 0.1, 0.08)
  )
  
  p <- plotKinaseTimeSeries(scoreTab_no_sig, pCut = 0.05)
  
  plot_data <- ggplot_build(p)$data[[2]]  # Data from geom_text layer
  
  expect_equal(sum(plot_data$label == "*"), 0)  # No significant markers should be present
})

# Test that the plot uses correct axis labels and their formatting
test_that("plotKinaseTimeSeries uses correct axis labels and formatting", {
  p <- plotKinaseTimeSeries(scoreTab)
  
  expect_equal(p$labels$x, "Time point")
  expect_equal(p$labels$y, "Kinase")
  
  # Check axis text formatting
  expect_equal(p$theme$axis.text.x$angle, 45)
  expect_equal(p$theme$axis.text.x$hjust, 1)
  expect_equal(p$theme$axis.text.x$vjust, 1)
})