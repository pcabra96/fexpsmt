################################################################################
##----------------------------------------------------------------------------##
## INDEX                                                                      ##
##----------------------------------------------------------------------------##
################################################################################

# 1. PACKAGES
# 2. SEED
# 3. SIMULATION PARAMETERS
# 4. SIMULATION
# 5. RESULTS
# 5.1. TIME
# 5.2. AVERAGE
# 5.3. TIME DOMAIN PARAMETER
# 5.4. FREQUENCY DOMAIN PARAMETER
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT

################################################################################
##----------------------------------------------------------------------------##
## 1. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

library(latex2exp)
library(gt)
library(magrittr)
library(htmltools)
library(ggplot2)
library(tidyr)
library(reshape2)

################################################################################
##----------------------------------------------------------------------------##
## 2. SEED                                                                    ##
##----------------------------------------------------------------------------##
################################################################################

set.seed(0)

################################################################################
##----------------------------------------------------------------------------##
## 3. SIMULATION PARAMETERS                                                   ##
##----------------------------------------------------------------------------##
################################################################################

PROCESS = "LMEMORY"
SUBPROCESS = "multiple"
symbol = "d"
d_coef_vec = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
# PATH TO LOAD THE DATA
################################################################################

path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/4. ",PROCESS,"/",SUBPROCESS,"/")

################################################################################
# 5.1 Running time
################################################################################

# DATA LISTS
own_times_farima_list = readRDS(file = paste0(path,"own_times_farima_list.RData"))
own_times_fexp_list = readRDS(file = paste0(path,"own_times_fexp_list.RData"))
r_times_list = readRDS(file = paste0(path,"r_times_list.RData"))

# DATA MATRICES
own_times_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_times_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_times_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(d_coef_vec)) {
  own_times_farima_matrix[wich_d,] = colMeans(own_times_farima_list[[wich_d]])
  own_times_fexp_matrix[wich_d,] = colMeans(own_times_fexp_list[[wich_d]])
  r_times_matrix[wich_d,] = colMeans(r_times_list[[wich_d]])
}

# fexpmst farima
rownames(own_times_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_times_farima_matrix) = paste0(POWER,"")
own_times_farima_matrix = as.data.frame(own_times_farima_matrix)

process_string = "FARIMA(0,d,0)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("fexpmst average running time (s) of ", N_SIMULATIONS, " {y<sub>", process_string, "<sub><sub>t</sub></sub></sub>}<sub>",equation_html)
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_times_farima_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 1.png",path = path)

# fexpmst fexp
rownames(own_times_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_times_fexp_matrix) = paste0(POWER,"")
own_times_fexp_matrix = as.data.frame(own_times_fexp_matrix)

process_string = "FEXP(0,d)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("fexpmst average running time (s) of ", N_SIMULATIONS, " {y<sub>", process_string, "<sub><sub>t</sub></sub></sub>}<sub>",equation_html)
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_times_fexp_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 2.png",path = path)

# fracdiff table
rownames(r_times_matrix) = paste0(d_coef_vec,"")
colnames(r_times_matrix) = paste0(POWER,"")
r_times_matrix = as.data.frame(r_times_matrix)

process_string = "FARIMA(0,d,0)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("fracdiff average running time (s) for ", N_SIMULATIONS, " {y<sub>", process_string, "<sub><sub>t</sub></sub></sub>}<sub>",equation_html)
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- r_times_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 3.png",path = path)

# COMPARISON 1 BIGGER
comparison_matrix_1 = r_times_matrix>own_times_farima_matrix
comparison_matrix_2 = r_times_matrix>own_times_fexp_matrix

comparison_matrix_1 = ifelse(comparison_matrix_1==TRUE,1,0)
df <- as.data.frame(comparison_matrix_1)
df$row_names <- rownames(comparison_matrix_1)
df_long <- pivot_longer(df, cols = -row_names, names_to = "col_names", values_to = "value")
col_order <- c("7", "8", "9", "10", "11", "12", "13", "14")
df_long$col_names <- factor(df_long$col_names, levels = col_order)

# Plot the matrix
ggplot(df_long, aes(x = col_names, y = row_names, fill = factor(value))) +
  geom_tile() +
  scale_fill_manual(values = c("white", "black"),labels = c("FALSE", "TRUE")) +
  labs(x = "T", y = "d", fill = "")+
  scale_x_discrete(labels = names)+
  ggtitle(TeX(paste0("Average running time: $fracdiff_{FARIMA}$ > $fexpmst_{FARIMA}$")))+
  theme(axis.title.y = element_text(angle = 0,vjust = 0.5), title = element_text(hjust = 0.5))

# Visualize the comparison
df <- as.data.frame(r_times_matrix/own_times_farima_matrix)
df$row_names <- rownames(comparison_matrix_1)
df_long <- pivot_longer(df, cols = -row_names, names_to = "col_names", values_to = "value")
col_order <- c("7", "8", "9", "10", "11", "12", "13", "14")
df_long$col_names <- factor(df_long$col_names, levels = col_order)

# Plot the matrix
ggplot(df_long, aes(x = col_names, y = row_names, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1,breaks = c(0,1, 5,10,15)) +
  labs(x = "T", y = "d")+
  scale_x_discrete(labels = names)+
  ggtitle(TeX(paste0("Average running time: $fexpmst_{FARIMA}$/$fracdiff_{FARIMA}$")))+
  theme(axis.title.y = element_text(angle = 0,vjust = 0.5), title = element_text(hjust = 0.5))

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("average running time (s): fracdiff<sub>FARIMA</sub>/fexpmst<sub>FARIMA</sub>")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (r_times_matrix/own_times_farima_matrix) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl

gtsave(gt_tbl, filename = "Table 4.png",path = path)

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("average running time (s): fracdiff<sub>FARIMA</sub>/fexpmst<sub>FEXP</sub>")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (r_times_matrix/own_times_fexp_matrix) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 5.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_times_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_times_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(own_times_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_times_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_times_farima_matrix = xtable(own_times_farima_matrix, digits = 5)
own_times_fexp_matrix = xtable(own_times_fexp_matrix, digits = 5)
r_times_matrix = xtable(r_times_matrix, digits = 5)

print(own_times_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(own_times_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(r_times_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

################################################################################
# 5.1 Running time
################################################################################

# DATA LISTS
own_average_farima_list = readRDS(file = paste0(path,"own_average_farima_list.RData"))
own_average_fexp_list = readRDS(file = paste0(path,"own_average_fexp_list.RData"))
r_average_list = readRDS(file = paste0(path,"r_average_list.RData"))

# DATA MATRICES
own_average_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_average_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_average_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(d_coef_vec)) {
  own_average_farima_matrix[wich_d,] = colMeans(own_average_farima_list[[wich_d]])
  own_average_fexp_matrix[wich_d,] = colMeans(own_average_fexp_list[[wich_d]])
  r_average_matrix[wich_d,] = colMeans(r_average_list[[wich_d]])
}

# fexpmst farima
rownames(own_average_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_average_farima_matrix) = paste0(POWER,"")
own_average_farima_matrix = as.data.frame(own_average_farima_matrix)

process_string = "FARIMA(0,d,0)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("fexpmst average fitted d ", N_SIMULATIONS, " {y<sub>", process_string, "<sub><sub>t</sub></sub></sub>}<sub>",equation_html)
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_average_farima_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl

#gtsave(gt_tbl, filename = "Table 1.png",path = path)

# fexpmst fexp
rownames(own_average_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_average_fexp_matrix) = paste0(POWER,"")
own_average_fexp_matrix = as.data.frame(own_average_fexp_matrix)

process_string = "FARIMA(0,d,0)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("fexpmst average fitted d ", N_SIMULATIONS, " {y<sub>", process_string, "<sub><sub>t</sub></sub></sub>}<sub>",equation_html)
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_average_fexp_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl

#gtsave(gt_tbl, filename = "Table 1.png",path = path)

# fracdiff fexp
rownames(r_average_matrix) = paste0(d_coef_vec,"")
colnames(r_average_matrix) = paste0(POWER,"")
r_average_matrix = as.data.frame(r_average_matrix)

process_string = "FARIMA(0,d,0)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("fexpmst average fitted d ", N_SIMULATIONS, " {y<sub>", process_string, "<sub><sub>t</sub></sub></sub>}<sub>",equation_html)
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- r_average_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl

#gtsave(gt_tbl, filename = "Table 1.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_average_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_average_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(own_average_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_average_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(r_average_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_average_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix_1 = as.matrix(abs(own_times_farima_matrix)<=abs(r_average_matrix))
comparison_matrix_2 = as.matrix(abs(own_average_fexp_matrix)<=abs(r_average_matrix))

################################################################################
# MSE COEFFICIENTS for d
################################################################################

# DATA
own_long_param_farima_list = readRDS(file = paste0(path,"own_long_param_farima_list.RData"))
own_long_param_fexp_list = readRDS(file = paste0(path,"own_long_param_fexp_list.RData"))
r_long_param_list = readRDS(file = paste0(path,"r_long_param_list.RData"))

# DATA MATRICES
own_long_param_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_long_param_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_long_param_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(d_coef_vec)) {
  own_long_param_farima_matrix[wich_d,] = colMeans( (own_long_param_farima_list[[wich_d]]-d_coef_vec[wich_d])^2 )
  own_long_param_fexp_matrix[wich_d,] = colMeans((own_long_param_fexp_list[[wich_d]]-d_coef_vec[wich_d])^2)
  r_long_param_matrix[wich_d,] = colMeans((r_long_param_list[[wich_d]]-d_coef_vec[wich_d])^2)
}




# AXIS NAMES FOR LATEX
rownames(own_long_param_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_long_param_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(own_long_param_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_long_param_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(r_long_param_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_long_param_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix_1 = as.matrix(own_long_param_farima_matrix<=r_long_param_matrix)
comparison_matrix_1 = as.matrix(own_long_param_fexp_matrix<=r_long_param_matrix)

# LATEX OUTPUT
own_long_param_farima_matrix = xtable(own_long_param_farima_matrix, digits = 5)
own_long_param_fexp_matrix = xtable(own_long_param_fexp_matrix, digits = 5)
r_long_param_matrix = xtable(r_long_param_matrix, digits = 5)

print(own_long_param_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(own_long_param_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(r_long_param_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x})

################################################################################
# MSE COEFFICIENTS for lambda=1
################################################################################

# DATA
own_lambda_farima_list = readRDS(file = paste0(path,"own_lambda_farima_list.RData"))
own_lambda_fexp_list = readRDS(file = paste0(path,"own_lambda_fexp_list.RData"))
r_lambda_list = readRDS(file = paste0(path,"r_lambda_list.RData"))

# DATA MATRICES
own_lambda_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_lambda_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_lambda_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(d_coef_vec)) {
  own_lambda_farima_matrix[wich_d,] = colMeans( (own_lambda_farima_list[[wich_d]]-1)^2 )
  own_lambda_fexp_matrix[wich_d,] = colMeans((own_lambda_fexp_list[[wich_d]]-1)^2)
  r_lambda_matrix[wich_d,] = colMeans((r_lambda_list[[wich_d]]-1)^2)
}

# AXIS NAMES FOR LATEX
rownames(own_lambda_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_lambda_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(own_lambda_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_lambda_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(r_lambda_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_lambda_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix_1 = as.matrix(own_lambda_farima_matrix<=r_lambda_matrix)
comparison_matrix_2 = as.matrix(own_lambda_fexp_matrix<=r_lambda_matrix)

# LATEX OUTPUT
own_lambda_farima_matrix = xtable(own_lambda_farima_matrix, digits = 5)
own_lambda_fexp_matrix = xtable(own_lambda_fexp_matrix, digits = 5)
r_lambda_matrix = xtable(r_lambda_matrix, digits = 5)

print(own_lambda_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(own_lambda_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(r_lambda_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x})

################################################################################
# Correct p-values for EXP(1)
################################################################################

# DATA
own_exp_1_farima_list = readRDS(file = paste0(path,"own_exp_1_farima_list.RData"))
own_exp_1_fexp_list = readRDS(file = paste0(path,"own_exp_1_fexp_list.RData"))
r_exp_1_list = readRDS(file = paste0(path,"r_exp_1_list.RData"))

# DATA MATRICES
own_exp_1_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_exp_1_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_exp_1_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (i in 1:length(d_coef_vec)) {
  own_exp_1_farima_matrix[i,] = colMeans(own_exp_1_farima_list[[i]]>=0.05)
  own_exp_1_fexp_matrix[i,] = colMeans(own_exp_1_fexp_list[[i]]>=0.05)
  r_exp_1_matrix[i,] = colMeans(r_exp_1_list[[i]]>=0.05)
}

# AXIS NAMES FOR LATEX
rownames(own_exp_1_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_exp_1_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(own_exp_1_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_exp_1_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

rownames(r_exp_1_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_exp_1_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

comparison_matrix_1 = as.matrix(own_exp_1_farima_matrix>=r_exp_1_matrix)
comparison_matrix_2 = as.matrix(own_exp_1_fexp_matrix>=r_exp_1_matrix)

# LATEX OUTPUT
own_exp_1_farima_matrix = xtable(own_exp_1_farima_matrix, digits = 5)
own_exp_1_fexp_matrix = xtable(own_exp_1_fexp_matrix, digits = 5)
r_exp_1_matrix = xtable(r_exp_1_matrix, digits = 5)

print(own_exp_1_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(own_exp_1_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)
print(r_exp_1_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_long_param_farima_matrix)), sanitize.text.function = function(x) {x})
