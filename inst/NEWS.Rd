\name{NEWS}
\title{News for Package \pkg{eRm}}
\encoding{UTF-8}
\section{Changes in Version 0.16-0}{
  \itemize{
	  \item \code{PersonMisfit()} incl. S3 added.
	}
}
\section{Changes in Version 0.15-8}{
  \itemize{
	 \item sanity check for missings in split vector in \code{LRtest.Rm()} added.
	 \item smoother confidence bands added in \code{plotGOF()}.
	}
}
\section{Changes in Version 0.15-7}{
  \itemize{
	  \item \code{LRtest.Rm()}: bugfix due to updates in base R.
	}
}
\section{Changes in Version 0.15-6}{
  \itemize{
	  \item \code{plotGOF()}: added arguments \code{x_axis}, \code{y_axis}, \code{set_par}, and \code{reset_par} to ease customization of the plot.
		\item Imports functions from default packages packages, as required by the new CRAN-check.
	}
}
\section{Changes in Version 0.15-5}{
  \itemize{
    \item implemented an \code{anova()} method for all models (except \acronym{LLRA}s, which have their own: \code{anova.llra}). see \code{?anova.eRm}
    \item added a function to compute separation reliability, as proposed in Wright & Stone (1999; see \code{?SepRel}).
    \item \code{plotINFO()} accepts more arguments via \code{\dots} (see \code{?plotINFO}).
    \item fixed a bug in \code{plotPWmap()}.
    \item fixed a bug in the internal function \code{get_item_cats()} related to \code{NA}s (affecting \acronym{LLRA} and item-information functions).
    \item switched encoding to \acronym{UTF}-8 to avoid problems with diacritics, etc. (e.g., Martin-Löf).
		\item updated citations.
    \item general improvements.
  }
}
\section{Changes in Version 0.15-4}{
  \itemize{
    \item streamlining some functions, updating the vignette, etc.
    \item warnings are now treated more consistently and can be suppressed with \code{suppressWarnings()}
    \item the \code{plotGOF()} function was revamped and tidied up in a couple of places, e.g.:
    \itemize{
      \item the plotting order of elements in a non-interactive call has been changed to put confidence ellipses and lines in the background and text in the foreground.
      \item the \eqn{x}{x}- and \eqn{y}{y}-limits are now dynamically computed by default, so that confidence regions and ellipses are inside the plotting region.
      \item the leading \dQuote{\code{I}} before item numbers has been removed for better legibility.
    }
    \item moved \code{NEWS} to the new fancy \code{NEWS.Rd} file/format
  }
}
\section{Changes in eRm version 0.15-3}{
  \itemize{
    \item an error in \code{test_info()} was fixed.
    \item \pkg{eRm} now depends on \code{R} \eqn{\geq3.0.0}{>= 3.0.0}
  }
}
\section{Changes in eRm version 0.15-2}{
  \itemize{
    \item an error in \code{Tpbis.stat} was fixed.
  }
}
\section{Changes in eRm version 0.15-1}{
  \itemize{
    \item a bug in \code{person.parameter()} has been fixed that caused the estimation to crash in some cases.
    \item a bug in \code{thresholds()} has been fixed that caused the routine to crash.
  }
}
\section{Changes in eRm version 0.15-0}{
  \itemize{
    \item the \pkg{RaschSampler} package has been merged into \pkg{eRm} for convenience (still available as a separate package).
  }
}
\section{Changes in eRm version 0.14-5}{
  \itemize{
    \item the package is now byte-compiled by default.
    \item some statistics added to \code{NPtest()}, \code{T7} and \code{T7a} removed.
    \item fixed a bug in \code{plotPWmap()}.
    \item fixed the \code{mplot} argument in plotting routines.
    \item fixed the split-criterion \code{"all.r"} in \code{LRtest()}.
    \item deleted all usages of \code{data()} in examples, since eRm uses \dQuote{lazy data.}
  }
}
\section{Changes in eRm version 0.14-4}{
  \itemize{
    \item when calling \code{NPtest()}, the \pkg{RaschSampler} can now be controlled more specifically (\code{burn_in}, \code{step}, \code{seed}).
    \item various improvements and bugfixes for LLRA-related functions.
    \item person parameter values can be extracted now for all persons using \code{coef()}.
      Additionally, in-/exclusion of extrapolated values (for 0 and perfect scores) can be controlled via the argument \code{extrapolated}.
    \item LRtest now computes standard errors (\code{se = TRUE}) by default.
    \item plotDIF now plots \dQuote{difficulties} for all models (formerly, \dQuote{easiness} parameters were plotted for Rasch models).
  }
}
\section{Changes in eRm version 0.14-3}{
  \itemize{
    \item minor bug fixed in \code{plotGOF()}, where on rare occasions confidence ellipses were plotted together with control lines (spotted by Peter Parker)
    \item improved labelling in \code{plotjointICC()}
  }
}
\section{Changes in eRm version 0.14-2}{
  \itemize{
    \item warning regarding group assignment when using median or mean split removed from \code{MLoef()}
    \item modification in \code{NPtest()} to split long output lines
    \item changed the delimiters of \code{plotDIF()} confidence intervals to \code{pch = 20} (small bullet).
  }
}
\section{Changes in eRm version 0.14-1}{
  \itemize{
    \item new experimental functions to calculate and plot item and test information (by Thomas Rusch)
    \item bug fixed in the calculation of item and person Infit t and Outfit t (hint from Rainer Alexandrowicz).
    \item \pkg{eRm} no longer depends on the \pkg{RaschSampler} package.
      However, it must be installed to use \code{NPtest()}.
    \item changed the delimiters of \code{plotDIF()} confidence intervals to \code{pch = 20}.
  }
}
\section{Changes in eRm version 0.14-0}{
  \itemize{
    \item new (wrapper) function \code{LLRA()} for fitting linear logistic models with relaxed assumptions including utilities for preparing data (\code{llra.datprep()}), setting up (\code{build_W()}) and modifying (\code{collapse_W()}) design matrices, comparing LLRA models (\code{anova()}) and plotting results (\code{plotTR()} and \code{plotGR()}) (by Thomas Rusch).
    \item \dQuote{exact} version of the Martin-Löf test for binary items and arbitrary splits added as method to \code{NPtest()}.
    \item in \code{plotGOF()} confidence ellipses can now be drawn for subsets of items, optionally using different colours
    \item new function \code{plotDIF()} (by Kathrin Gruber): plots confidence intervals for item parameters estimated separately in subgroups, uses LR objects as input
    \item adapted the \code{MLoef()} function to work with polytomous data and more than two item groups
    \item error checks in NPtest:
    \enumerate{
      \item 0/full responses for items meaningless for NPtest,
      \item group in \code{method = "T4"} must be of type logical,
      \item specifying all items for T4 gives meaningless results.
    }
    \item warning regarding group assignment when using median split removed from \code{LRtest()} and \code{Waldtest()}.
    \item some modifications in \code{plotPWmap()}: horizontal plotting, different default plotting symbols, option to change size of plotting symbols
    \item bug in \code{MLoef()} fixed (now using logs in calculating the person contributions)
    \item \pkg{eRm} now depends on \code{R} \eqn{\geq2.12.0}{>= 2.12.0}
    \item Latin1 encoding removed
    \item bug in \code{plotICC()} (always same title) fixed
  }
}
\section{Changes in eRm version 0.13-0}{
  \itemize{
    \item \code{LLTM()}, \code{LRSM()}, and \code{LPCM()} work now for repeated measurement designs with treatment groups and missing values.
    \item Rename vignette to \dQuote{eRm}.
  }
}
\section{Changes in eRm version 0.12-2}{
  \itemize{
    \item new function \code{plotPWmap()} to plot Bond-and-Fox style pathway maps for the data by Julian Gilbey.
      Since calculation of the \eqn{t}{t}-statistics requires calculation of the kurtosis of the standardized residuals, according changes to \code{itemfit.ppar()}, \code{personfit.ppar()}, \code{pifit.internal()}, \code{print.ifit()}, and \code{print.pfit()}.
    \item \code{plotPImap()} patched by Julian Gilbey: length of \code{item.subset} did not match the documentation, warning stars did not all appear, pre-calculated person.parameter data can be passed to the function via pp, mis-ordered items can be coloured.
      some minor bugs fixed.
    \item the optimizer can be changed to \code{optim()} using \code{fitctrl <- "optim"} and reset to \code{nlm()} (the default) with \code{fitctrl <- "nlm"}
    \item value of \code{LRtest()} now contains the list \code{fitobj} which contains the model objects according to the subgroups specified by \code{splitcr}
    \item \code{MLoef()} no longer supports missing values
  }
}
\section{Changes in eRm version 0.12-1}{
  \itemize{
    \item function invalid from package \pkg{gtools} integrated into \pkg{eRm}.
      \pkg{eRm} no longer depends on \pkg{gtools}.
  }
}
\section{Changes in eRm version 0.12-0}{
  \itemize{
    \item for \code{RM()}, \code{RSM()}, and \code{PCM()}: eta parameters are now displayed as difficulty parameters; \code{print()} and \code{summary()} methods changed accordingly.
    \item new labeling of eta parameters in \code{RM()}, \code{RSM()}, and \code{PCM()}.
      they now are labeled according to the estimated parameters for items (\code{RM()}), items + categories (\code{RSM()}), items x categories (\code{PCM()})
    \item function \code{MLoef()} for Martin-Löf-Test added
    \item \code{df} in \code{personfit()} and \code{itemfit()} corrected
    \item the \code{logLik()} functions now extract the log-likelhood and df into objects of class \code{"logLik.eRm"} and \code{"loglik.ppar"} with elements \code{loglik} and \code{df}.
      the corresponding print methods have been modified accordingly.
    \item method \code{coef.ppar()} to extract person parameter estimates added
    \item option for beta parameters added to coef.eRm
    \item in confint.eRm: default \code{parm = "beta"}
    \item minor modifications in the help file for \code{IC()}
    \item \code{plotPImap()}: revised rug added, bug concerning \code{item.subset} fixed, minor modifications to enhance readability
    \item minor modifications in \code{plotjointICC()}: allows for main title and colors, option \code{legpos = FALSE} suppresses legends, dev.new removed, \code{legend = FALSE} produced incorrect labeling
    \item minor modifications in \code{plotICC()}: allows for main title and colors, default coloring with \code{col = NULL} instead of \code{NA} for compatibility, option \code{legpos = FALSE} suppresses legends, \code{mplot} is now \code{FALSE} if only one item is specified
    \item plot.ppar: dev.new removed
    \item option \dQuote{visible} in print.ifit und print.pfit to allow for avoiding overly long output and for extraction of infit and outfit values (maybe changed to a coef method later)
    \item \code{strwrap()} for NPtest print methods to break long lines
    \item new methods \code{IC.default()} and \code{pmat.default()} for enhanced error messages
    \item lazy loading package and datafiles
  }
}
