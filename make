\documentclass[twoside,a4paper]{report}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{amsmath, amssymb, graphics, setspace , bbm}
\usepackage{listings}
\usepackage{color}
\usepackage{verbatim}
\usepackage{float}
\usepackage{xstring}
\usepackage{hyperref}

%\usepackage{polski}
%\usepackage[utf8]{inputenc}
%\usepackage[OT4]{fontenc}

\usepackage{stylePHD}

% Global style for lstlisting .f and .f90 code.
% When entering .f code use:
% \lstset{language = [77]Fortran}
% \begin{lstlisting}
% ...
% \end{lstlisting}
% When using .f90 code use:
% \lstset{language = [90]Fortran}
% \begin{lstlisting}
% ...
% \end{lstlisting}

\lstset{
	backgroundcolor=\color[rgb]{0.9 , 0.9 , 0.9},
	basicstyle=\footnotesize\ttfamily,
	keywordstyle=\bfseries\color[rgb]{0.5 , 0.5 , 0.5},
	commentstyle=\itshape\color[rgb]{1,0 , 0.0 , 0.0},
	identifierstyle=\color[rgb]{0.0 , 0.0 , 0.0},	
	stringstyle=\color[rgb]{1.0 , 0.0 , 1.0}
}


\title{The two-nucleon and three-nucleon systems in three dimensions.}
\author{Kacper Topolnicki}

\begin{document}

\maketitle

\tableofcontents

\input{abstract.tex}

\input{Preface/Preface.tex}

\input{notation/notation.tex}

\input{introduction/introduction.tex}

\input{degrees_basic/degrees_basic.tex}

\input{scattering/scattering.tex}

\input{bound_state_fad/bound_state_fad.tex}

\input{degrees/degrees.tex}

\input{detail_2N/detail_2N.tex}

\input{transitionoperator/transitionoperator.tex}

\input{currents/currents.tex}

\input{detail_3N/detail_3N.tex}

\input{numerical_methods/numerical_methods.tex}

\input{util1N2N3N/util1N2N3N.tex}

\input{FunctionArray/FunctionArray.tex}

\input{CodeOrganization/CodeOrganization.tex}

\input{appendix/appendix.tex}

\input{link_to_PWD/link_to_PWD.tex}

\bibliographystyle{unsrt}
\bibliography{bibl}

\end{document}
