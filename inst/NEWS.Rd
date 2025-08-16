% borrowed some formatting from
% http://svn.r-project.org/R/trunk/doc/NEWS.Rd
% -*- coding: utf-8 -*-
%%\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}
%%\newcommand{\Rlogo}{\if{html}{\figure{../../html/logo.jpg}}\if{latex}{\figure{logo.jpg}{options:
%width=0.5in}}}
\newcommand{\GITpkg}{\href{https://github.com/abarbour/#1}{github}}
\newcommand{\JK}{\href{https://github.com/jkennel}{@jkennel}}
\newcommand{\CRANpkg}{\href{https://cran.r-project.org/package=#1}{CRAN}}
\newcommand{\kitagawa}{\pkg{kitagawa}}
\newcommand{\gitkit}{\GITpkg{kitagawa} page}

\name{kitagawa}
\title{
News for package \kitagawa{}: \CRANpkg{kitagawa} -- \GITpkg{kitagawa}
}

\encoding{UTF-8}

\section{CHANGES IN VERSION 3.1+ (last updated 2025-08-16)}{
  \itemize{
   
    \item 3.1.3: added link to psd to satisfy CRAN error 
    
    \item CI is handled via Github Actions (replaces Travis-CI)
    
    \item 3.1.2: Fixed error in vignette engine specification
    
    \item 3.1.1: added a \code{ticks.only} option to \code{logticks}

    \item 3.1: Added the model of Wang et al (2018; \href{https://doi.org/10.1029/2018WR022793}{10.1029/2018WR022793}),
    to the model toolbox and to the ResponseModels toolbox. Thanks to \JK for contributing this.
    \item New test validating Wang model against Hsieh model (they are equivalent when specific leakage is zero).
    Thanks to \JK for contributing this.
    \item Previous build on CRAN broke because sapa/ifultools were removed; hence, the specifics of
    \code{cross_spectrum} have changed. The code now uses new functionality from the \CRANpkg{psd} package. 
    (sapa/ifultools were unreliablely maintained for a while.)
    One benefit is this new approach allows for adaptive refinement of the (sine-multitaper) cross spectrum rather
    than uniform tapering across all frequencies. But, one unfortunate consequence is
    that any codes/scripts using the Welch option (i.e., \code{k=NULL}) will break.
  }
}

\section{CHANGES IN VERSION 2.3 (2019-09-18)}{
  \itemize{
    \item \href{https://github.com/abarbour/kitagawa/pull/2/commits/1610331a03a4470fa83e30160af4a6c279072277}{Fixes} a sign error in the Liu solution. Thanks to \JK for catching this.
    \item Removed the need for pandoc-citeproc in Tohoku vignette because of a presently unfixed bug.
  }
}

\section{CHANGES IN VERSION 2.2-2 (2018-09-14)}{
  \itemize{
    \item NA values stop execution of \code{cross_spectrum}; and, added tests.
    \item NEW FUNCTION: \code{logsmoo} which can be used to smooth hdr signals in log-space
    \item UPDATED: Tohoku cross-spectrum vignette
    \item Added DOIs for references
  }
}

\section{CHANGES IN VERSION 2.2-1 (2018-01-20)}{
  \itemize{
    \item Added optional switch to \code{cross_spectrum} to switch from multitaper to Welch overlapping estimation
  }
}

\section{CHANGES IN VERSION 2.2-0 (2017-07-11)}{
  \itemize{
    \item Tohoku cross-spectrum vignette
    \item Function for cross-spectrum computation
    \item Includes Windows CI via AppVeyor
    \item unit-testing framework (testthat, \CRANpkg{testthat}), with some tests
  }
}


\section{CHANGES IN VERSION 2.1-1 (2015-02-10)}{
  \itemize{
    \item Mostly a maintenance release.
    \item Travis-CI is now used to track builds; see also the \gitkit{}.
  }
}

\section{CHANGES IN VERSION 2.1-0 (2013-10-03)}{
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

\section{CHANGES IN VERSION 1.1-1 (2013-03-12)}{
  \itemize{
    \item Removed empty \code{.Last.lib} function.
    \item Added missing closures in 'NEWS.Rd'.
  }
}

\section{CHANGES IN VERSION 1.1-0 (2013-03-11)}{
  \itemize{
    \item Added \code{.onAttach} message into 'zzz.R' file.
    \item NEWS is now in Rd format in the inst/ directory.
  }
}

\section{CHANGES IN VERSION 1.0-1 (2012-10-25)}{
  \itemize{
    \item This marks the first release to \CRANpkg{kitagawa}.
    \item added a demo
    \item bug fix: when given a single frequency (omega), the calculation failed
    \item documentation improved
  }
}

\section{CHANGES IN VERSION 1.0-0 (2012-09-27)}{
  \itemize{
    \item First version, checks ok.
    \item Documentation complete, thanks to roxygen2 (\CRANpkg{roxygen2}).
  }
}
