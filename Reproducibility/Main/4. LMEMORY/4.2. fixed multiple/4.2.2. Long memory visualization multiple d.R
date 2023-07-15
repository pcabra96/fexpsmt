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
# 5.1. TIME (Tables: 1,2,3,4,5)
# 5.2. AVERAGE (Tables: 6,7,8,9,10)
# 5.3. TIME DOMAIN PARAMETER d (Tables: 11,12,13,14,15)
# 5.4. FREQUENCY DOMAIN PARAMETER (Tables: 16,17,18,19,20)
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT (Tables: 21,22,23,24,25)

################################################################################
##----------------------------------------------------------------------------##
## 1. PACKAGES                                                                ##
##----------------------------------------------------------------------------##
################################################################################

library(RColorBrewer)
library(latex2exp)
library(gt)
library(gtExtras)
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
SUBPROCESS = "fixed multiple"
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/4. ",PROCESS,"/4.2. ",SUBPROCESS,"/")
symbol = "d"
d_coef_vec = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))

################################################################################
# 5.1. TIME (Tables: 1,2,3,4,5)
################################################################################

# LOAD DATA LISTS
own_times_farima_list = readRDS(file = paste0(path,"own_times_farima_list.RData"))
own_times_fexp_list = readRDS(file = paste0(path,"own_times_fexp_list.RData"))
r_times_list = readRDS(file = paste0(path,"r_times_list.RData"))

# COMPUTE AVERAGES: DATA MATRICES
own_times_farima_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
own_times_fexp_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))
r_times_matrix = matrix(0,nrow = length(d_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(d_coef_vec)) {
  own_times_farima_matrix[wich_d,] = colMeans(own_times_farima_list[[wich_d]])
  own_times_fexp_matrix[wich_d,] = colMeans(own_times_fexp_list[[wich_d]])
  r_times_matrix[wich_d,] = colMeans(r_times_list[[wich_d]])
}

#-------------------------------------------------------------------------------
# 5.1.1. fexpmst farima
#-------------------------------------------------------------------------------

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

# AXIS NAMES FOR LATEX
rownames(own_times_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_times_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_times_farima_matrix = xtable(own_times_farima_matrix, digits = 5)
print(own_times_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)

#Back to normal
rownames(own_times_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_times_farima_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.1.2. fexpmst fexp
#-------------------------------------------------------------------------------

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

# AXIS NAMES FOR LATEX
rownames(own_times_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_times_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_times_fexp_matrix = xtable(own_times_fexp_matrix, digits = 5)
print(own_times_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x},booktabs = TRUE)

#Back to normal
rownames(own_times_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_times_fexp_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.1.3. fracdiff farima
#-------------------------------------------------------------------------------

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

# AXIS NAMES FOR LATEX
rownames(r_times_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_times_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
r_times_matrix = xtable(r_times_matrix, digits = 5)
print(r_times_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

#Back to normal
rownames(r_times_matrix) = paste0(d_coef_vec,"")
colnames(r_times_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.1.4. fracdiff farima vs fexpmst farima
#-------------------------------------------------------------------------------

pal <- function(x){
  f_neg <- scales::col_numeric(palette = c(rev(c("white",brewer.pal(9, "OrRd")[1:3]))),domain = c(0,1))
  f_pos <- scales::col_numeric(palette = c("white",brewer.pal(8, "Blues")[3:5]), domain = c(1, 17))
  ifelse(x < 1, f_neg(x), f_pos(x))}

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("average running time (s): fracdiff<sub>FARIMA</sub>/fexpmst<sub>FARIMA</sub>")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (r_times_matrix/own_times_farima_matrix) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |>
          fmt_number(columns = everything(), decimals = 5) |>
          data_color(columns = everything(), colors = pal)
gt_tbl

gtsave(gt_tbl, filename = "Table 4.png",path = path)

#-------------------------------------------------------------------------------
# 5.1.5. fracdiff farima vs fexpmst fexp
#-------------------------------------------------------------------------------

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("average running time (s): fracdiff<sub>FARIMA</sub>/fexpmst<sub>FEXP</sub>")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (r_times_matrix/own_times_fexp_matrix) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)|>
          data_color(columns = everything(), colors = pal)

gt_tbl
gtsave(gt_tbl, filename = "Table 5.png",path = path)

################################################################################
# 5.2. AVERAGE (Tables: 6,7,8,9,10)
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
  own_average_farima_matrix[wich_d,] = colMeans((own_average_farima_list[[wich_d]])^2)
  own_average_fexp_matrix[wich_d,] = colMeans((own_average_fexp_list[[wich_d]])^2)
  r_average_matrix[wich_d,] = colMeans((r_average_list[[wich_d]])^2)
}

#-------------------------------------------------------------------------------
# 5.2.1. fexpmst farima
#-------------------------------------------------------------------------------

rownames(own_average_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_average_farima_matrix) = paste0(POWER,"")
own_average_farima_matrix = as.data.frame(own_average_farima_matrix)

process_string = "FARIMA(0,d,0)"
package_string = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(y&#x0304<sub>fracdiff<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_average_farima_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 6.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_average_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_average_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_average_farima_matrix = xtable(own_average_farima_matrix, digits = 5)
print(own_average_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(own_average_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_average_farima_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.2.2. fexpmst fexpmst
#-------------------------------------------------------------------------------

rownames(own_average_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_average_fexp_matrix) = paste0(POWER,"")
own_average_fexp_matrix = as.data.frame(own_average_fexp_matrix)

process_string = "FEXP(0,d)"
package_string = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(y&#x0304<sub>fracdiff<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_average_fexp_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 7.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_average_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_average_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_average_fexp_matrix = xtable(own_average_fexp_matrix, digits = 5)
print(own_average_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(own_average_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_average_fexp_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.2.3. fracdiff farima
#-------------------------------------------------------------------------------

rownames(r_average_matrix) = paste0(d_coef_vec,"")
colnames(r_average_matrix) = paste0(POWER,"")
r_average_matrix = as.data.frame(r_average_matrix)

process_string = "FARIMAP(0,d,0)"
package_string = "fracdiff"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(y&#x0304<sub>fracdiff<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- r_average_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 8.png",path = path)

# AXIS NAMES FOR LATEX
rownames(r_average_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_average_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
r_average_matrix = xtable(r_average_matrix, digits = 5)
print(r_average_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(r_average_matrix) = paste0(d_coef_vec,"")
colnames(r_average_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.2.4. fexmpst farima vs fracdiff farima
#-------------------------------------------------------------------------------

pal <- function(x){
  f_neg <- scales::col_numeric(palette = c(rev(c("white",brewer.pal(8, "OrRd")[1:3]))),domain = c(-2,1))
  f_pos <- scales::col_numeric(palette = c("white",brewer.pal(9, "Blues")[2:5]), domain = c(1, 9))
  ifelse(x < 1, f_neg(x), f_pos(x))}

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(y&#x0304<sub>fracdiff<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)/MSE(y&#x0304<sub>fexpmst<sub><sub>FARIMAP(0,d,0)</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (abs(r_average_matrix)/abs(own_average_farima_matrix)) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5) |>
          data_color(columns = everything(), colors = pal)

gt_tbl

gtsave(gt_tbl, filename = "Table 9.png",path = path)

# AXIS NAMES FOR LATEX
rownames(r_average_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_average_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
r_average_matrix = xtable(r_average_matrix, digits = 5)
print(r_average_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(r_average_matrix) = paste0(d_coef_vec,"")
colnames(r_average_matrix) = paste0(POWER,"")

comparison_matrix_1 = as.matrix(abs(own_times_farima_matrix)<=abs(r_average_matrix))
comparison_matrix_2 = as.matrix(abs(own_average_fexp_matrix)<=abs(r_average_matrix))

#-------------------------------------------------------------------------------
# 5.2.5. fexmpst farima vs fracdiff farima
#-------------------------------------------------------------------------------

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(y&#x0304<sub>fracdiff<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)/MSE(y&#x0304<sub>fexpmst<sub><sub>FEXP(0,d)</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (abs(r_average_matrix)/abs(own_average_fexp_matrix)) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)|>
          data_color(columns = everything(), colors = pal)

gt_tbl
gtsave(gt_tbl, filename = "Table 10.png",path = path)

# AXIS NAMES FOR LATEX
rownames(r_average_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_average_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
r_average_matrix = xtable(r_average_matrix, digits = 5)
print(r_average_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(r_average_matrix) = paste0(d_coef_vec,"")
colnames(r_average_matrix) = paste0(POWER,"")

################################################################################
# 5.3. TIME DOMAIN PARAMETER d (Tables: 11,12,13,14,15)
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

#-------------------------------------------------------------------------------
# 5.3.1. fexpmst farima
#-------------------------------------------------------------------------------

rownames(own_long_param_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_long_param_farima_matrix) = paste0(POWER,"")
own_long_param_farima_matrix = as.data.frame(own_long_param_farima_matrix)

process_string = "FARIMA(0,d,0)"
package_string = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted d<sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_long_param_farima_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 11.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_long_param_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_long_param_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_long_param_farima_matrix = xtable(own_long_param_farima_matrix, digits = 5)
print(own_long_param_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(own_long_param_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_long_param_farima_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.3.2. fexpmst fexpmst
#-------------------------------------------------------------------------------

rownames(own_long_param_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_long_param_fexp_matrix) = paste0(POWER,"")
own_long_param_fexp_matrix = as.data.frame(own_long_param_fexp_matrix)

process_string = "FEXP(0,d)"
package_string = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted d<sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_long_param_fexp_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 12.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_long_param_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_long_param_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_long_param_fexp_matrix = xtable(own_long_param_fexp_matrix, digits = 5)
print(own_long_param_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(own_long_param_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_long_param_fexp_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.3.3. fexpmst fexpm
#-------------------------------------------------------------------------------

rownames(r_long_param_matrix) = paste0(d_coef_vec,"")
colnames(r_long_param_matrix) = paste0(POWER,"")
r_long_param_matrix = as.data.frame(r_long_param_matrix)

process_string = "FEXP(0,d)"
package_string = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted d<sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- r_long_param_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 13.png",path = path)

# AXIS NAMES FOR LATEX
rownames(r_long_param_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_long_param_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
r_long_param_matrix = xtable(r_long_param_matrix, digits = 5)
print(r_long_param_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(r_long_param_matrix) = paste0(d_coef_vec,"")
colnames(r_long_param_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.3.4. fexmpst farima vs fracdiff farima
#-------------------------------------------------------------------------------

pal <- function(x){
  f_neg <- scales::col_numeric(palette = c(rev(c("white",brewer.pal(10, "OrRd")[2:4]))),domain = c(0.85,1))
  f_pos <- scales::col_numeric(palette = c("white",brewer.pal(10, "Blues")[2:4]), domain = c(1, 1.3))
  ifelse(x < 1, f_neg(x), f_pos(x))}

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted d<sub>fracdiff<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)/MSE(fitted d<sub>fexpmst<sub><sub>FARIMAP(0,d,0)</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (abs(r_long_param_matrix)/abs(own_long_param_farima_matrix)) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5) |>
          data_color(columns = everything(), colors = pal)

gt_tbl
gtsave(gt_tbl, filename = "Table 14.png",path = path)

#-------------------------------------------------------------------------------
# 5.3.5. fexmpst farima vs fracdiff farima
#-------------------------------------------------------------------------------

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted d<sub>fracdiff<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)/MSE(fitted d<sub>fexpmst<sub><sub>FEXP(0,d)</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (abs(r_long_param_matrix)/abs(own_long_param_fexp_matrix)) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5) |>
          data_color(columns = everything(), colors = pal)

gt_tbl
gtsave(gt_tbl, filename = "Table 15.png",path = path)

################################################################################
# 5.4. FREQUENCY DOMAIN PARAMETER (Tables: 16,17,18,19,20)
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

#-------------------------------------------------------------------------------
# 5.4.1. fexpmst farima
#-------------------------------------------------------------------------------

rownames(own_lambda_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_lambda_farima_matrix) = paste0(POWER,"")
own_lambda_farima_matrix = as.data.frame(own_lambda_farima_matrix)

process_string = "FARIMA(0,d,0)"
package_string = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted &lambda;<sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_lambda_farima_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 16.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_lambda_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_lambda_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_lambda_farima_matrix = xtable(own_long_param_farima_matrix, digits = 5)
print(own_lambda_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(own_lambda_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_lambda_farima_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.4.2. fexpmst fexp
#-------------------------------------------------------------------------------

rownames(own_lambda_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_lambda_fexp_matrix) = paste0(POWER,"")
own_lambda_fexp_matrix = as.data.frame(own_lambda_fexp_matrix)

process_string = "FEXP(0,d)"
package_string = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted &lambda;<sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_lambda_fexp_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 17.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_lambda_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_lambda_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_lambda_fexp_matrix = xtable(own_long_param_farima_matrix, digits = 5)
print(own_lambda_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(own_lambda_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_lambda_fexp_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.4.3. fracdiff farima
#-------------------------------------------------------------------------------

rownames(r_lambda_matrix) = paste0(d_coef_vec,"")
colnames(r_lambda_matrix) = paste0(POWER,"")
r_lambda_matrix = as.data.frame(r_lambda_matrix)

process_string = "FARIMA(0,d,0)"
package_string = "fracdiff"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted &lambda;<sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- r_lambda_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 18.png",path = path)

# AXIS NAMES FOR LATEX
rownames(r_lambda_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_lambda_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
r_lambda_matrix = xtable(own_long_param_farima_matrix, digits = 5)
print(r_lambda_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(r_lambda_matrix) = paste0(d_coef_vec,"")
colnames(r_lambda_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.4.4. fexmpst farima vs fracdiff farima
#-------------------------------------------------------------------------------

pal <- function(x){
  f_neg <- scales::col_numeric(palette = c(rev(c("white",brewer.pal(5, "OrRd")[1:3]))),domain = c(0.7,1))
  f_pos <- scales::col_numeric(palette = c("white",brewer.pal(3, "Blues")[1:3]), domain = c(1, 1.3))
  ifelse(x < 1, f_neg(x), f_pos(x))}

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted &lambda;<sub>fracdiff<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)/MSE(fitted &lambda;<sub>fexpmst<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (abs(r_lambda_matrix)/abs(own_lambda_farima_matrix)) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5) |>
          data_color(columns = everything(), colors = pal)

gt_tbl
gtsave(gt_tbl, filename = "Table 19.png",path = path)

#-------------------------------------------------------------------------------
# 5.4.5. fexmpst farima vs fracdiff farima
#-------------------------------------------------------------------------------

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted &lambda;<sub>fracdiff<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)/MSE(fitted &lambda;<sub>fexpmst<sub><sub>FEXP(0,d)</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (abs(r_lambda_matrix)/abs(own_lambda_fexp_matrix)) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)|>
          data_color(columns = everything(), colors = pal)

gt_tbl
gtsave(gt_tbl, filename = "Table 20.png",path = path)

################################################################################
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT (Tables: 21,22,23,24,25)
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

#-------------------------------------------------------------------------------
# 5.5.1. fexpmst farima
#-------------------------------------------------------------------------------

rownames(own_exp_1_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_exp_1_farima_matrix) = paste0(POWER,"")
own_exp_1_farima_matrix = as.data.frame(own_exp_1_farima_matrix)

process_string = "FARIMA(0,d,0)"
package_string = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("Percentage of non rejected H<sub>0</sub>: I<sup>*</sup><sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub> &sim; exp(&lambda;=1)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_exp_1_farima_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 21.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_exp_1_farima_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_exp_1_farima_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_exp_1_farima_matrix = xtable(own_long_param_farima_matrix, digits = 5)
print(own_exp_1_farima_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(own_exp_1_farima_matrix) = paste0(d_coef_vec,"")
colnames(own_exp_1_farima_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.5.2. fexpmst fexp
#-------------------------------------------------------------------------------

rownames(own_exp_1_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_exp_1_fexp_matrix) = paste0(POWER,"")
own_exp_1_fexp_matrix = as.data.frame(own_exp_1_fexp_matrix)

process_string = "FEXP(0,d)"
package_string = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("Percentage of non rejected H<sub>0</sub>: I<sup>*</sup><sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub> &sim; exp(&lambda;=1)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- own_exp_1_fexp_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 22.png",path = path)

# AXIS NAMES FOR LATEX
rownames(own_exp_1_fexp_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(own_exp_1_fexp_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
own_exp_1_fexp_matrix = xtable(own_long_param_farima_matrix, digits = 5)
print(own_exp_1_fexp_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(own_exp_1_fexp_matrix) = paste0(d_coef_vec,"")
colnames(own_exp_1_fexp_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.5.3. fexpmst fexp
#-------------------------------------------------------------------------------

rownames(r_exp_1_matrix) = paste0(d_coef_vec,"")
colnames(r_exp_1_matrix) = paste0(POWER,"")
r_exp_1_matrix = as.data.frame(r_exp_1_matrix)

process_string = "FARIMA(0,d,0)"
package_string = "fracdiff"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("Percentage of non rejected H<sub>0</sub>: I<sup>*</sup><sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub> &sim; exp(&lambda;=1)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- r_exp_1_matrix |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_times_matrix)) |>
  cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
  tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)

gt_tbl
gtsave(gt_tbl, filename = "Table 23.png",path = path)

# AXIS NAMES FOR LATEX
rownames(r_exp_1_matrix) = paste0("$\\mathbf{",d_coef_vec,"}$")
colnames(r_exp_1_matrix) = paste0("$\\mathbf{2^{",POWER,"}}$")

# LATEX OUTPUT
r_exp_1_matrix = xtable(own_long_param_farima_matrix, digits = 5)
print(r_exp_1_matrix, include.rownames = TRUE, hline.after = c(-1,0, nrow(own_times_farima_matrix)), sanitize.text.function = function(x) {x})

# Back to normal
rownames(r_exp_1_matrix) = paste0(d_coef_vec,"")
colnames(r_exp_1_matrix) = paste0(POWER,"")

#-------------------------------------------------------------------------------
# 5.5.4. fexmpst farima vs fracdiff farima
#-------------------------------------------------------------------------------

pal <- function(x){
  f_neg <- scales::col_numeric(palette = c(rev(c("white",brewer.pal(3, "OrRd")[1:3]))),domain = c(0, 0.8,1))
  f_pos <- scales::col_numeric(palette = c("white",brewer.pal(3, "Blues")[1:3]), domain = c(1, 1.2))
  ifelse(x < 1, f_neg(x), f_pos(x))}

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("(H<sub>0</sub>:<sub>fracdiff<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)/(H<sub>0</sub>:<sub>fexpmst<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (abs(r_exp_1_matrix)/abs(own_exp_1_farima_matrix)) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_exp_1_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)|>
          data_color(columns = everything(), colors = pal)

gt_tbl
gtsave(gt_tbl, filename = "Table 24.png",path = path)

#-------------------------------------------------------------------------------
# 5.5.5. fexmpst farima vs fracdiff farima
#-------------------------------------------------------------------------------

process_string = ""
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("(H<sub>0</sub>:<sub>fracdiff<sub><sub>FARIMA(0,d,0)</sub></sub></sub>)/(H<sub>0</sub>:<sub>fexpmst<sub><sub>FEXP(0,d)</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- (abs(r_exp_1_matrix)/abs(own_exp_1_fexp_matrix)) |> gt(rowname_col = "d", rownames_to_stub = TRUE) |> tab_header(title =  html(main)) |> tab_spanner(label = "T", columns = colnames(r_exp_1_matrix)) |>
          cols_label("7" = html("2<sup>7</sup>"),"8" = html("2<sup>8</sup>"),"9" = html("2<sup>9</sup>"),"10" = html("2<sup>10</sup>"),"11" = html("2<sup>11</sup>"),"12" = html("2<sup>12</sup>"),"13" = html("2<sup>13</sup>"),"14" = html("2<sup>14</sup>")) |>
          tab_stubhead(label = "d")|> opt_css(css) |> fmt_number(columns = everything(), decimals = 5)|>
          data_color(columns = everything(), colors = pal)

gt_tbl

gtsave(gt_tbl, filename = "Table 25.png",path = path,expand = 30)


