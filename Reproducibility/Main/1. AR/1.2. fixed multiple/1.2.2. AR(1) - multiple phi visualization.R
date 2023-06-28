################################################################################
# PACKAGES
################################################################################

library(latex2exp)
library(gt)
library(magrittr)
library(htmltools)
library(ggplot2)
library(tidyr)
library(reshape2)

################################################################################
# PARAMETERS
################################################################################

PROCESS = "AR"
SUBPROCESS = "fixed multiple"
ar_coef_vec = c(-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9)
POWER = 7:14
N_SIMULATIONS = 1000
names=c(TeX("$2^7$"), TeX("$2^8$"), TeX("$2^9$"), TeX("$2^{10}$"), TeX("$2^{11}$"), TeX("$2^{12}$"), TeX("$2^{13}$"),TeX("$2^{14}$"))
path = paste0("~/Documents/2. UNIGE/2023-1 Master Thesis/fexpsmt/Reproducibility/Main/1. ",PROCESS,"/1.2. ",SUBPROCESS,"/")

################################################################################
# 5.1. TIME
################################################################################

time_own_list = readRDS(file = paste0(path,"time_own_list.RData"))
time_r_list = readRDS(file = paste0(path,"time_r_list.RData"))

average_time_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
average_time_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  average_time_own[i,] = colMeans((time_own_list[[i]]))
  average_time_r[i,] = colMeans((time_r_list[[i]]))
}

#-------------------------------------------------------------------------------
# 5.1.1. fexpmst AR
#-------------------------------------------------------------------------------

param = average_time_r
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

package = "stats"
process_string = "AR(1)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0(package," average running time (s) of ", N_SIMULATIONS, " {y<sub>", process_string, "<sub><sub>t</sub></sub></sub>}<sub>",equation_html)
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 1.png",path = path)

#-------------------------------------------------------------------------------
# 5.1.2. stats AR
#-------------------------------------------------------------------------------

param = average_time_r
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

package = "stats"
process_string = "AR(1)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0(package," average running time (s) of ", N_SIMULATIONS, " {y<sub>", process_string, "<sub><sub>t</sub></sub></sub>}<sub>",equation_html)
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 2.png",path = path)


################################################################################
# 5.2. AVERAGE
################################################################################

# DATA LISTS
average_own_list = readRDS(file = paste0(path,"average_own_list.RData"))
average_r_list = readRDS(file = paste0(path,"average_r_list.RData"))

# DATA MATRICES
own_average_farima_matrix = matrix(0,nrow = length(ar_coef_vec), ncol = length(POWER))
r_average_matrix = matrix(0,nrow = length(ar_coef_vec), ncol = length(POWER))

for (wich_d in 1:length(ar_coef_vec)) {
  own_average_farima_matrix[wich_d,] = colMeans((average_own_list[[wich_d]])^2)
  r_average_matrix[wich_d,] = colMeans((average_r_list[[wich_d]])^2)
}

#-------------------------------------------------------------------------------
# 5.2.1. stats AR
#-------------------------------------------------------------------------------

param = own_average_farima_matrix
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

package = "fexpmst"
process_string = "AR(1)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(y&#x0304<sub>",package,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 3.png",path = path)

#-------------------------------------------------------------------------------
# 5.2.2 stats AR
#-------------------------------------------------------------------------------

param = r_average_matrix
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

package = "stats"
process_string = "AR(1)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(y&#x0304<sub>",package,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 4.png",path = path)

#-------------------------------------------------------------------------------
# 5.2.3. fexmpst AR vs stats AR
#-------------------------------------------------------------------------------

param = (r_average_matrix/own_average_farima_matrix)
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

process_string = "AR(1)"
package_1 = "stats"
package_2 = "fexpmst"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(y&#x0304<sub>",package_1,"<sub><sub>",process_string,"</sub></sub></sub>)/MSE(y&#x0304<sub>",package_2,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 5.png",path = path)


################################################################################
# 5.3. MSE COEFFICIENTS for phi
################################################################################

# DATA
fit_own_coef_list = readRDS(file = paste0(path,"fit_own_coef_list.RData"))
fit_r_coef_list = readRDS(file = paste0(path,"fit_r_coef_list.RData"))

# DATA MATRICES
MSE_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
MSE_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))

for (i in 1:length(ar_coef_vec)) {
  MSE_own[i,] = colMeans((fit_own_coef_list[[i]]-ar_coef_vec[i])^2)
  MSE_r[i,] = colMeans((fit_r_coef_list[[i]]-ar_coef_vec[i])^2)
}

#-------------------------------------------------------------------------------
# 5.3.1 fexpmst AR
#-------------------------------------------------------------------------------

param = MSE_own
package = "fexpmst"
main = paste0("MSE(fitted &varphi;<sub>",package,"<sub><sub>",process_string,"</sub></sub></sub>)")
process_string = "AR(1)"

rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 6.png",path = path)

#-------------------------------------------------------------------------------
# 5.3.2 stats AR
#-------------------------------------------------------------------------------

param = MSE_r
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

package = "stats"
process_string = "AR(1)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted &varphi;<sub>",package,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 7.png",path = path)

#-------------------------------------------------------------------------------
# 5.3.3 fexpmst vs stats AR
#-------------------------------------------------------------------------------

param = (MSE_r/MSE_own)
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

package_1 = "stats"
package_2 = "fexpmst"
process_string = "AR(1)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted &varphi;<sub>",package_1,"<sub><sub>",process_string,"</sub></sub></sub>)/MSE(fitted &varphi;<sub>",package_2,"<sub><sub>",package_2,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 8.png",path = path)

################################################################################
# 5.4. FREQUENCY DOMAIN PARAMETER
################################################################################

# DATA
fit_own_exp_list = readRDS(file = paste0(path,"fit_own_exp_list.RData"))
fit_r_exp_list = readRDS(file = paste0(path,"fit_r_exp_list.RData"))

# OWN CODE
MSE_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
MSE_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  MSE_own[i,] = colMeans((fit_own_exp_list[[i]]-1)^2)
  MSE_r[i,] = colMeans((fit_r_exp_list[[i]]-1)^2)
}

#-------------------------------------------------------------------------------
# 5.4.1 fexpmst AR
#-------------------------------------------------------------------------------

param = MSE_own
package_string = "fexpmst"
main = paste0("MSE(fitted &lambda;<sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub>)")
process_string = "AR(1)"

rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 9.png",path = path)

#-------------------------------------------------------------------------------
# 5.3.2 stats AR
#-------------------------------------------------------------------------------

param = MSE_r
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

package_string = "stats"
process_string = "AR(1)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted &lambda;<sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 10.png",path = path)

#-------------------------------------------------------------------------------
# 5.4.3 fexpmst vs stats AR
#-------------------------------------------------------------------------------

param = (MSE_r/MSE_own)
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

package_1 = "stats"
package_2 = "fexpmst"
process_string = "AR(1)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("MSE(fitted &lambda;<sub>",package_1,"<sub><sub>",process_string,"</sub></sub></sub>)/MSE(fitted &lambda;<sub>",package_2,"<sub><sub>",package_2,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 11.png",path = path)

################################################################################
# 5.5. FREQUENCY DOMAIN GOODNESS OF FIT
################################################################################

# DATA
p_val_own_exp_list = readRDS(file = paste0(path,"p_val_own_exp_list.RData"))
p_val_r_exp_list = readRDS(file = paste0(path,"p_val_r_exp_list.RData"))

# OWN CODE
prop_non_rejection_own = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
prop_non_rejection_r = matrix(0, nrow = length(ar_coef_vec), ncol = length(POWER))
for (i in 1:length(ar_coef_vec)) {
  prop_non_rejection_own[i,] = colMeans(p_val_own_exp_list[[i]]>=0.05)
  prop_non_rejection_r[i,] = colMeans(p_val_r_exp_list[[i]]>=0.05)
}

#-------------------------------------------------------------------------------
# 5.4.1 fexpmst AR
#-------------------------------------------------------------------------------

param = prop_non_rejection_own
package_string = "fexpmst"
process_string = "AR(1)"
main = paste0("Percentage of non rejected H<sub>0</sub>: I<sup>*</sup><sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub> &sim; exp(&lambda;=1)")

rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 12.png",path = path)

#-------------------------------------------------------------------------------
# 5.3.2 stats AR
#-------------------------------------------------------------------------------

param = prop_non_rejection_r
package_string = "stats"
process_string = "AR(1)"
main = paste0("Percentage of non rejected H<sub>0</sub>: I<sup>*</sup><sub>",package_string,"<sub><sub>",process_string,"</sub></sub></sub> &sim; exp(&lambda;=1)")

rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 13.png",path = path)

#-------------------------------------------------------------------------------
# 5.4.3 fexpmst vs stats AR
#-------------------------------------------------------------------------------

param = (prop_non_rejection_r/prop_non_rejection_own)
rownames(param) = paste0(ar_coef_vec,"")
colnames(param) = paste0(POWER,"")
param = as.data.frame(param)

package_1 = "stats"
package_2 = "fexpmst"
process_string = "AR(1)"
equation <- "T<br>t=1"
equation_html <- paste0("<div class='equation'>", equation, "</div>")
main = paste0("(H<sub>0</sub>:<sub>",package_1,"<sub><sub>",process_string,"</sub></sub></sub>)/(H<sub>0</sub>:<sub>",package_2,"<sub><sub>",process_string,"</sub></sub></sub>)")
css <- ".equation {display: inline-block;vertical-align: middle;}"

gt_tbl <- param |> gt(rowname_col = html("Phi"), rownames_to_stub = TRUE) |> tab_spanner(label = "T", columns = colnames(average_time_own) ) |>
  fmt_number(columns = everything(), decimals = 5) |>
  tab_stubhead(label = HTML("&varphi;<sub>1</sub>") ) |> tab_header(title =  HTML(main)) |> opt_css(css)|>
  cols_label("7" = HTML("2<sup>7</sup>"),"8" = HTML("2<sup>8</sup>"),"9" = HTML("2<sup>9</sup>"),"10" = HTML("2<sup>10</sup>"),"11" = HTML("2<sup>11</sup>"),"12" = HTML("2<sup>12</sup>"),"13" = HTML("2<sup>13</sup>"),"14" = HTML("2<sup>14</sup>"))

gt_tbl

gtsave(gt_tbl, filename = "Table 14.png",path = path)
