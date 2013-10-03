% borrowed some formatting from
% http://svn.r-project.org/R/trunk/doc/NEWS.Rd
% -*- coding: utf-8 -*-
%%\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}
%%\newcommand{\Rlogo}{\if{html}{\figure{../../html/logo.jpg}}\if{latex}{\figure{logo.jpg}{options:
%width=0.5in}}}
\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{kitagawa}
\title{News for package 'kitagawa'}
\encoding{UTF-8}

\section{CHANGES IN kitagawa VERSION 2.1-0 (2013-10-03)}{
  \itemize{
    \item New capability to model open well response: \code{open_well_response}
    \item \code{well_response} now returns complex frequency response only, rather than
    returnin amplitude and phase as separate vectors.
    \item Substantial documentation improvements (and shuffling around).
    \item Response functions now add attributes and class to output: 'wrsp' and 'owrsp'
    for sealed and open well response, respectively.
    \item New methods for 'wrsp' and 'owrsp', including \code{plot.wrsp}, for example.
    \item \code{kitplot} is simply a link to \code{plot.wrsp}
    \item Changed variable names in \code{well_response} and added defaults in case some
    are missing.
    \item Added a vignette highlighting model differences.
    \item New utility function \code{.in0to1} to test for existence between zero to one;
    this is useful to check if a parameter is valid (e.g. Skemptons ratio).
    \item New utility \code{constants} to access stored physical constants.
  }
}
\section{CHANGES IN kitagawa VERSION 1.1-1 (2013-03-12)}{
  \itemize{
    \item Removed empty \code{.Last.lib} function.
    \item Added missing closures in 'NEWS.Rd'.
  }
}
\section{CHANGES IN kitagawa VERSION 1.1-0 (2013-03-11)}{
  \itemize{
    \item Added \code{.onAttach} message into 'zzz.R' file.
    \item NEWS is now in Rd format in the inst/ directory.
  }
}
\section{CHANGES IN kitagawa VERSION 1.0-1 (2012-10-25)}{
  \itemize{
    \item This marks the first release to CRAN: \CRANpkg{kitagawa}.
    \item added a demo
    \item bug fix: when given a single frequency (omega), the calculation failed
    \item documentation improved
  }
}
\section{CHANGES in kitagawa VERSION 1.0-0 (2012-09-27)}{
  \itemize{
    \item First version, checks ok.
    \item Documentation complete, thanks to \CRANpkg{roxygen2}.
  }
}
