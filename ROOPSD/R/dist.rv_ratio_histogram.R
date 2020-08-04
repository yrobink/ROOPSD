
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

#' rv_ratio_histogram 
#'
#' @description
#' rv_ratio_histogram distribution in OOP way.
#'
#' @details
#' Fit separatly P( X < x | X > 0 ) and P(X=0)
#'
#' @examples
#' ## Generate sample
#' X = numeric(10000)
#' X[1:2000] = 0
#' X[2001:10000] = stats::rexp( n = 8000 , rate = 1 )
#'
#' ## And fit it
#' rvX = rv_ratio_histogram$new()
#' rvX$fit( X , x0 = 0 )
#'
#' @export
rv_ratio_histogram = R6::R6Class( "rv_ratio_histogram" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field rvXp [ROOPSD::rv_histogram] Describes P(X < x | X > x0)
	rvXp = NULL,
	#' @field x0 [double] location of mass: P( X = x0 )
	x0   = NULL,
	#' @field p0 [double] p0 = P( X = x0 )
	p0   = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new rv_ratio_histogram object.
    #' @param ... If a param `Y` and `x0` is given, the fit method is called with `...`.
    #' @return A new `rv_ratio_histogram` object.
	initialize = function(...)
	{
		kwargs = list(...)
		if( !is.null(kwargs[["Y"]]) && !is.null(kwargs[["x0"]]) )
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
	
	## cdf ##{{{
	#' @description
    #' Cumulative Distribution Function
    #' @param q [vector] Quantiles to compute the CDF
    #' @return cdf values
	cdf = function( q )
	{
		cdf = base::rep( 0 , length(q) )
		idxp = q > 0
		idx0 = !idxp
		cdf[idxp] = ( 1 - self$p0 ) * self$rvXp$cdf(q[idxp]) + self$p0
		cdf[idx0] = self$p0 / 2
		return(cdf)
	},
	##}}}
	
	## icdf ##{{{
	#' @description
    #' Inverse of Cumulative Distribution Function
    #' @param p [vector] Probabilities to compute the CDF
    #' @return icdf values
	icdf = function( p )
	{
		idxp = p > self$p0
		idx0 = !idxp
		icdf = base::rep( 0 , length(p) )
		icdf[idxp] = self$rvXp$icdf( (p[idxp] - self$p0) / ( 1 - self$p0 ) )
		icdf[idx0] = 0
		return(icdf)
	},
	##}}}
	
	## sf ##{{{
	#' @description
    #' Survival Function
    #' @param q [vector] Quantiles to compute the SF
    #' @return sf values
	sf = function( q ) 
	{
		return( 1. - self$cdf(q) )
	},
	##}}}
	
	## isf ##{{{
	#' @description
    #' Inverse of Survival Function
    #' @param p [vector] Probabilities to compute the SF
    #' @return isf values
	isf = function( p ) 
	{
		return(self$icdf(1. - p))
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit method for the histograms
    #' @param Y [vector] Dataset to infer the histogram
    #' @param x0 [double] Location of mass point
    #' @param bins [vector or integer] bins values
    #' @return `self`
	fit = function( Y , x0 , bins = as.integer(100) )
	{
		self$x0   = x0
		Yp = Y[Y>self$x0]
		self$rvXp = ROOPSD::rv_histogram$new( Y = Yp , bins = bins )
		self$p0 = base::sum(!(Y>x0)) / length(Y)
		
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
	
	
	
	)
)

