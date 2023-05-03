# fexpsmt
This will be an R package to simulate FARIMA and FEXP processes from their spectral density. It is part of the master thesis of the University of Geneva. The project is under the supervision of professor Davide La Vecchia.

{\Huge 1. Objective and definitions} \\

$\newline$
The goal of this master thesis is to create and implement an algorithm to simulate any univariate autoregressive fractionally integrated moving average (FARIMA) or fractionally exponential (FEXP) time series $\{Y_t\}_{t=1}^T$ from their theoretic spectral densities $f_Y(\omega)$. 

$\newline$
A process $\{Y_t\}_{t=1}^T$ is said to be FARIMA$(p,d,q)$ (also called ARFIMA) if it is stationary and it satisfies the difference equation \ref{eq:1} or equivalenty if it has the spectral density defined in equation \ref{eq:2}
(See \hyperlink{Beran}{[9] Beran, Jan. (1994)}).

\begin{equation}
\label{eq:1}
\phi(L)(1-L)^dY_t=\theta(L)\epsilon_t
\end{equation}

In equation \ref{eq:1} $L$ is the backshift operator, $\varepsilon_t$ is a gaussian noise with mean 0 and variance $\sigma^2$, $\theta(L)$ and $\phi(L)$ are polynomials of degree $p$ and $q$ respectively, and $d$ regulates the long memory component of the process. If $d=0$ the process has only short term behaviour, if $0<d<0.5$ the process has long memory and if $-0.5<d<0$ the process is antipersistent.

\begin{equation}
\label{eq:2}
f_{FARIMA}(\omega)= \frac{\sigma^2}{2\pi}|1-e^{i\omega}|^{-2d}\frac{|1-\sum_{j=1}^{q}\theta_je^{i\omega}|^2}{|1-\sum_{j=1}^{p}\phi_je^{i\omega}|^2}
\end{equation}

In equation \ref{eq:2} the parameter $d$ also regulates the long term behaviour, $\sigma^2$ is the variance of $\varepsilon_t$ of equation \ref{eq:1} and the order of the polynomials $p$ and $q$ is the same as the one of \ref{eq:1}.

$\newline$
This model is a generalization of the ARMA model, which was popularized by
\hyperlink{BoxJenkins}{[2] Box, G. E., Jenkins, G. M., $\And$ Reinsel, G. (1970)}.
It includes a long memory component which was first discussed in the context of hidrology 
(See \hyperlink{Hurst}{[1] Hurst, H. E. (1951)}). Some years after that it was formalized most notably by polish matematician Benoit Mandelbrot. It was finally included as an extension to the ARMA model thanks to the research of 
\hyperlink{Granger}{[4] Granger, C. W., $\And$ Joyeux, R. (1980)}.

$\newline$
A process $\{Y_t\}_{t=1}^T$ is said to be FEXP if it has the spectral density defined in equation \ref{eq:3} (See \hyperlink{Beran}{[9] Beran, Jan. (1994)}).

\begin{equation}
\label{eq:3}
f_{FEXP}(\omega)= \frac{\sigma^2}{2\pi}|1-e^{i\omega}|^{-2d}exp\{2\sum_{k=0}^{K}\eta_kcos(k\omega)\}
\end{equation}

$\newline$
This model is a generalization of the exponential model (EXP) proposed by 
\hyperlink{Bloomfield}{[3] Bloomfield, P. (1973)}. It comes from the observation that the logarithm of the estimated spectral density $log(\hat{f}(\omega))$ is a well behaved function that can be represented as a truncated fourier series.

$\newline$
{\Huge 2. Motivation} \\

$\newline$
This thesis is motivated by three main factors. Firstly, there is a deep interest in spectral analysis in time series. Secondly, there is a desire to increase computing efficiency. While it is possible to design an algorithm to simulate the FEXP model in the time domain, it is believed that using the spectral domain in conjunction with the Fast Fourier Transform, which has a complexity of $\mathcal{O}(Tlog(T))$, can lead to lower complexity. Thirdly, there is a possibility to build a package in R that could be published in CRAN.

$\newline$
To achieve the objectives of this thesis, two models were chosen for study. The FARIMA process was selected because it is a popular parametric model that has been widely used by academics and practitioners since its introduction. Furthermore, there are already good packages in R to simulate such processes, making it possible to access a benchmark to measure the quality of the algorithm.

$\newline$
The FEXP model was chosen because a simulation of such a model has not been developed yet in any known statistical software, making it a novelty. Additionally, this model is closely related to the FARIMA model. In fact, a FEXP model with infinite parameters recovers a FARIMA process exactly, which means that it is also possible to check the correctness of the algorithm by comparing the results of the two models.

$\newline$
Overall, the combination of these models and techniques will provide a valuable contribution to the field of time series analysis, and the potential to publish a package in CRAN provides an added incentive to ensure the quality and accuracy of the research.

$\newpage$

{\Huge 3. Theoretical framework} \\

$\newline$
In order to check the correctness of the simulation method for the FEXP process we will use the fact that the log-FEXP is the truncated Fourier Series of log-spectral density of the FARIMA process. Hence, if the method is correct for the FEXP process, it should converge to the FARIMA process as parameters go to infinity. Let us prove that relationship.

Let us start with two alternative definitions of the ARMA process with equation \ref{eq:4} and equation \ref{eq:5}.

\begin{equation}
\label{eq:4}
log(f_{ARMA}(\omega))= log(\frac{\sigma^2}{2\pi})+\sum_{k=1}^{\infty}[\sum_{j=1}^{p}\frac{b_j^k}{k}-\sum_{j=1}^{q}\frac{a_j^k}{k}]2cos(k\omega)
\end{equation}

Where $b_j$ and $a_j$ are roots of the polynomials $1-\sum_{j=1}^{p}\phi_je^{i\omega}$ and $1-\sum_{j=1}^{q}\theta_je^{i\omega}$ respectively.

\begin{equation}
\label{eq:5}
log(f_{ARMA}(\omega))= log(\frac{\sigma^2}{2\pi})+2\sum_{k=1}^{\infty}c_kcos(k\omega)
\end{equation}

Note that $log(f_{ARMA}(\omega) = log(f_{ARMA}(\omega+2\pi))$. As every periodic function it can be expressed as a fourier series, as in equation \ref{eq:6}:

\begin{equation}
\label{eq:6}
log(f_{ARMA}(\omega)) = A_0 + \sum_{n=1}^{\infty}A_ncos(n\omega)+\sum_{n=1}^{\infty}B_nsin(n\omega)
\end{equation}

Where:
$\newline$
- $A_0 = \frac{1}{2\pi} \int_{-\pi}^{\pi}log(f_{ARMA}(\omega))d\omega$
$\newline$
- $A_n = \frac{1}{\pi} \int_{-\pi}^{\pi}log(f_{ARMA}(\omega))cos(n\omega)d\omega$
$\newline$
- $B_n = \frac{1}{\pi} \int_{-\pi}^{\pi}log(f_{ARMA}(\omega))sin(n\omega)d\omega$

$\newline$
We also know the following properties:
$\newline$
- $sin(t) = -sin(t)$, so it is an odd function.
$\newline$
- $cos(t) = cos(-t)$, so it is an even function.
$\newline$
- $log(f_{ARMA}(\omega)=log(f_{ARMA}(-\omega))$, so it is an even function.

$\newline$
Therefore, we can conclude that $B_n=0 \ \forall n=[1,2,...\infty]$ because the symmetric integral of an odd function is zero. Hence equation \ref{eq:6} can be rewritten as equation \ref{eq:7}:

\begin{equation}
\label{eq:7}
log(f_{ARMA}(\omega)) = A_0 + \sum_{n=1}^{\infty}A_ncos(n\omega)
\end{equation}

Since the Fourier Series representation of a function is unique due to Cantor's uniqueness theorem (See \hyperlink{CantorTheorem}{[6] J. Marshall Ash (1989)}), one can conclude the following about equations \ref{eq:3}, \ref{eq:4}, \ref{eq:5} and \ref{eq:7}:
$\newline$
- $A_0=log(\frac{\sigma^2}{2\pi})$
$\newline$
- $A_n = 2[\sum_{j=1}^{p}\frac{b_j^k}{k}-\sum_{j=1}^{q}\frac{a_j^k}{k}] = 2c_k = 2\eta_k \ \forall n,k \in \mathbb{N} \ | \ n=k$

This implies that the truncated Fourier Series representation of the log-spectral density of an ARMA process is corresponds to the EXP of order $p$. It is trivial to show that the same can be stated about the FARIMA and the FEXP.
$\newline$

{\Huge 4. Algorithm} \\

This algorithm is in the spirit of \hyperlink{Percival}{[8] Percival, D. B. (1993)}. 
It has been used also by other researchers in multivariate setting (See \hyperlink{Chambers}{[11]Chambers, M. J. (1995)}) and functional time-series analysis 
(See \hyperlink{Rubin}{[14] Rub\'{i}n, T., $\And$ Panaretos, V. M. (2020)}). 
It has some minor changes and is general in the sense that It can be used whenever one has some theoretical spectral density of stationary processes.

\begin{algorithm}
\caption{Time Series Simulation}
\begin{algorithmic}[1]

\State $\textbf{Size of the process}$
\State $M \ \gets \ 2\times T$

\State $\textbf{Theoretical spectral densit at specific frequencies}$
\State $\omega_t \ \gets \  \frac{2\pi}{M} t, \ \ t=0,1,2...,\frac{M}{2}$

\State $\textbf{Simulate Gaussian Noise}$
\State $\mathcal{W}_t \  \gets \  \mathcal{N}(0,1), \ \ t=0,1,2...,M-1$

\State $\textbf{Induce randomness to the theoretical spectral density}$
\State $\mathcal{V}_0 \  \gets \  \sqrt{f_Y(\omega_0)}\mathcal{W}_0$

\State $\mathcal{V}_t \  \gets \  \sqrt{\frac{1}{2}f_Y(\omega_t)}(\mathcal{W}_t+i\mathcal{W}_{2t-1}), \ \ 1 \leq t< \frac{M}{2}$

\State $\mathcal{V}_{\frac{M}{2}} \  \gets \  \sqrt{f_Y(\omega_{\frac{M}{2}})}\mathcal{W}_{\frac{M}{2}}$

\State $\mathcal{V}_t \  \gets \  \mathcal{V^*}_{M-t}, \ \ \frac{M}{2} < t \leq M-1$ \Comment{the asterisk denotes complex conjugation}

\State $\textbf{From frequency to time domain}$
\State $Y_{1:T} \  \gets \   \frac{1}{\sqrt{M}} iFFT(\mathcal{V}_{0:M-1})$ \Comment{iFFT is the inverse Fast Fourier Transform}

\end{algorithmic}
\end{algorithm}

$\newline$

{\Huge 5. Initial Results} \\

To this point some succesful simulations of FARIMA process have been performed. Since one know the true parameters, the success is measured with the mean squared error of the estimated FARIMA parameters against the true parameters.

$\newpage$

{\Huge References} \\

$\newline$
- \hypertarget{Hurst}{[1] Hurst, H. E. (1951)}. 
Long-term storage capacity of reservoirs. Transactions of the American society of civil engineers, 116(1), 770-799.
$\newline$
- \hypertarget{BoxJenkins}{[2] Box, G. E., Jenkins, G. M., $\And$ Reinsel, G. (1970)}. Time series analysis: forecasting and control Holden-day San Francisco. BoxTime Series Analysis: Forecasting and Control Holden Day1970.
$\newline$
- \hypertarget{Bloomfield}{[3] Bloomfield, P. (1973)}. An exponential model for the spectrum of a scalar time series. Biometrika, 60(2), 217-226.
$\newline$
- \hypertarget{Granger}{[4] Granger, C. W., $\And$ Joyeux, R. (1980)}. An introduction to long‐memory time series models and fractional differencing. Journal of time series analysis, 1(1), 15-29.
$\newline$
- \hypertarget{Davies}{[5] Davies, R. B., $\And$ Harte, D. S. (1987)}. Tests for Hurst effect. Biometrika, 74(1), 95-101.
$\newline$
- \hypertarget{CantorTheorem}{[6] J. Marshall Ash (1989)}. Uniqueness of Representation by Trigonometric Series, The American Mathematical Monthly, 96:10, 873-885, DOI: 10.1080/00029890.1989.11972299
$\newline$
- \hypertarget{Brockwell}{[7] Brockwell, P. J., & Davis, R. A. (1991)}. Time series: theory and methods, 2nd ed., Springer-Verlag, New York.
$\newline$
- \hypertarget{Percival}{[8] Percival, D. B. (1993)}. Simulating Gaussian random processes with specified spectra. Computing Science and Statistics, 534-534.
$\newline$
- \hypertarget{Beran}{[9] Beran, Jan. (1994)}. Statistics for Long-Memory Processes. Vol. 61. CRC press.
$\newline$
- \hypertarget{Hamilton}{[10] Hamilton, J. D. (1994)}. Time Series Analysis. Princeton University Press. 
$\newline$
- \hypertarget{Chambers}{[11]Chambers, M. J. (1995)}. The simulation of random vector time series with given spectrum. Mathematical and computer modelling, 22(2), 1-6.
$\newline$
- \hypertarget{ssssss}{[12] Shumway, R. H., Stoffer, D. S., $\And$ Stoffer, D. S. (2000)}. Time series analysis and its applications (Vol. 3). New York: springer.
$\newline$
- \hypertarget{History}{[13] Graves, T., Gramacy, R., Watkins, N., $\And$ Franzke, C. (2017)}. A brief history of long memory: Hurst, Mandelbrot and the road to ARFIMA, 1951–1980. Entropy, 19(9), 437.
$\newline$
- \hypertarget{Rubin}{[14] Rub\'{i}n, T., $\And$ Panaretos, V. M. (2020)}. Spectral simulation of functional time series. arXiv preprint arXiv:2007.08458.
