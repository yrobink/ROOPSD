
################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2020                                                ##
##                                                                            ##
## yoann.robin.k@gmail.com                                                    ##
##                                                                            ##
## This software is a computer program that is part of the ROOPSD (R Object   ##
## Oriented Programming for Statistical Distribution). This library makes it  ##
## possible to draw, fit and access to CDF and SF function                    ##
## for classic statistical distribution.                                      ##
##                                                                            ##
## This software is governed by the CeCILL-C license under French law and     ##
## abiding by the rules of distribution of free software.  You can  use,      ##
## modify and/ or redistribute the software under the terms of the CeCILL-C   ##
## license as circulated by CEA, CNRS and INRIA at the following URL          ##
## "http://www.cecill.info".                                                  ##
##                                                                            ##
## As a counterpart to the access to the source code and  rights to copy,     ##
## modify and redistribute granted by the license, users are provided only    ##
## with a limited warranty  and the software's author,  the holder of the     ##
## economic rights,  and the successive licensors  have only  limited         ##
## liability.                                                                 ##
##                                                                            ##
## In this respect, the user's attention is drawn to the risks associated     ##
## with loading,  using,  modifying and/or developing or reproducing the      ##
## software by the user in light of its specific status of free software,     ##
## that may mean  that it is complicated to manipulate,  and  that  also      ##
## therefore means  that it is reserved for developers  and  experienced      ##
## professionals having in-depth computer knowledge. Users are therefore      ##
## encouraged to load and test the software's suitability as regards their    ##
## requirements in conditions enabling the security of their systems and/or   ##
## data to be ensured and,  more generally, to use and operate it in the      ##
## same conditions as regards security.                                       ##
##                                                                            ##
## The fact that you are presently reading this means that you have had       ##
## knowledge of the CeCILL-C license and that you accept its terms.           ##
##                                                                            ##
################################################################################
################################################################################

################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2020                                                ##
##                                                                            ##
## yoann.robin.k@gmail.com                                                    ##
##                                                                            ##
## Ce logiciel est un programme informatique faisant partie de la librairie   ##
## ROOPSD (R Object Oriented Programming for Statistical Distribution). Cette ##
## librarie permet de tirer aléatoirement, inférer et acceder aux fonctions   ##
## CDF et SF des distributions statistisques classiques.                      ##
##                                                                            ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et  ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez    ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions     ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA   ##
## sur le site "http://www.cecill.info".                                      ##
##                                                                            ##
## En contrepartie de l'accessibilité au code source et des droits de copie,  ##
## de modification et de redistribution accordés par cette licence, il n'est  ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le     ##
## titulaire des droits patrimoniaux et les concédants successifs.            ##
##                                                                            ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques      ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au      ##
## développement et à la reproduction du logiciel par l'utilisateur étant     ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à      ##
## manipuler et qui le réserve donc à des développeurs et des professionnels  ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les    ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du     ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la       ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,   ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         ##
##                                                                            ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez     ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les  ##
## termes.                                                                    ##
##                                                                            ##
################################################################################
################################################################################

## rv_histogram ##{{{

#' rv_histogram 
#'
#' @description
#' rv_histogram distribution in OOP way.
#'
#' @details
#' Use quantile to fit the histogram
#'
#' @examples
#' ## Generate sample
#' X = numeric(10000)
#' X[1:5000] = stats::rnorm( n = 5000 , mean = 2 , sd = 1 )
#' X[5000:10000] = stats::rexp( n = 5000 , rate = 1 )
#'
#' ## And fit it
#' rvX = rv_histogram$new()
#' rvX$fit(X)
#'
#' @export
rv_histogram = R6::R6Class( "rv_histogram" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field min [double] min value for the estimation
	min   = NULL,
	#' @field max [double] max value for the estimation
	max   = NULL,
	#' @field tol [double] numerical tolerance
	tol = 1e-2,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new rv_histogram object.
    #' @param ... If a param `Y` is given, the fit method is called with `...`.
    #' @return A new `rv_histogram` object.
	initialize = function(...)
	{
		kwargs = list(...)
		if( !is.null(kwargs[["Y"]]) )
			base::do.call( self$fit , kwargs )
	},
	##}}}
	
	## rvs ##{{{
	#' @description
    #' Generation sample from the histogram
    #' @param n [integer] Number of samples drawn
    #' @return A vector of samples
	rvs = function( n )
	{
		p = stats::runif( n , 0 , 1 )
		return(self$icdf(p))
	},
	##}}}
	
	## density ##{{{
	#' @description
    #' Density function
    #' @param x [vector] Values to compute the density
    #' @return density
	density = function( x ) 
	{
		return(private$.density(x))
	},
	##}}}
	
	## logdensity ##{{{
	#' @description
    #' Log density function
    #' @param x [vector] Values to compute the log-density
    #' @return the log density
	logdensity = function( x ) 
	{
		return(base::log(private$.density(x)))
	},
	##}}}
	
	## cdf ##{{{
	#' @description
    #' Cumulative Distribution Function
    #' @param q [vector] Quantiles to compute the CDF
    #' @return cdf values
	cdf = function( q )
	{
		return(private$.cdf(q))
	},
	##}}}
	
	## icdf ##{{{
	#' @description
    #' Inverse of Cumulative Distribution Function
    #' @param p [vector] Probabilities to compute the CDF
    #' @return icdf values
	icdf = function( p )
	{
		return(private$.icdf(p))
	},
	##}}}
	
	## sf ##{{{
	#' @description
    #' Survival Function
    #' @param q [vector] Quantiles to compute the SF
    #' @return sf values
	sf = function( q ) 
	{
		return( 1. - private$.cdf(q) )
	},
	##}}}
	
	## isf ##{{{
	#' @description
    #' Inverse of Survival Function
    #' @param p [vector] Probabilities to compute the SF
    #' @return isf values
	isf = function( p ) 
	{
		return(private$.icdf(1. - p))
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit method for the histograms
    #' @param Y [vector] Dataset to infer the histogram
    #' @param bins [vector or integer] bins values
    #' @return `self`
	fit = function( Y , bins = as.integer(1000) )
	{
		self$min = base::min(Y)
		self$max = base::max(Y)
		
		## Start with cdf and icdf
		Yr = base::rank( Y , ties.method = "max" )
		Ys = base::sort(Y)
		Yru = base::sort(base::unique(Yr))
		p  = Yru / length(Y)
		q  = Ys[Yru]
		
		p[1] = 0
#		p = base::c( 0 , p )
#		q = base::c( self$min , q )
#		if( q[1] == q[2] )
#		{
#			eps  = base::sqrt(.Machine$double.xmin)
#			q[2] = (1-eps) * q[1] + eps * q[3]
#		}
#		q = base::c( self$min - .Machine$double.xmin , q )
		
		private$.cdf  = stats::approxfun( q , p , method = "linear" , ties = "ordered" , yleft = 0 , yright = 1 )
		private$.icdf = stats::approxfun( p , q , method = "linear" , ties = "ordered" , yleft = NaN , yright = NaN )
		
		## Continue with density
		x = NULL
		if( is.integer(bins) )
		{
			delta = 1e-2 * (self$max - self$min)
			x = base::seq( self$min - delta , self$max + delta , length = bins )
		}
		else
		{
			x = bins
		}
		
		## Density function
		hist = graphics::hist( Y , breaks = x , plot = FALSE )
		p = hist$density #/ base::sum(hist$density)
		c = hist$mids
		private$.density = stats::approxfun( c , p ,
		                                      yleft = 0 ,
		                                      yright = 0 ,
		                                      ties = "ordered" )
		
		return(self)
	}
	##}}}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	.cdf     = NULL,
	.icdf    = NULL,
	.density = NULL
	
	
	)
)
##}}}


