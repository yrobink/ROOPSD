
################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2021                                                ##
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
## Copyright Yoann Robin, 2021                                                ##
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

#' mrv_histogram 
#'
#' @description
#' Multivariate rv_histogram distribution in OOP way.
#'
#' @details
#' Used for a multivariate dataset, fit each marge
#'
#' @examples
#' ## Generate sample
#' X = matrix( stats::rnorm( n = 10000 ) , ncol = 4 )
#'
#' ## And fit it
#' rvX = mrv_histogram$new()
#' rvX$fit(X)
#'
#' @export
mrv_histogram = R6::R6Class( "mrv_histogram" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field n_features [integer] Number of features (dimensions)
	n_features = NULL,
	#' @field law_ [list] List of marginal distributions
	law_       = NULL,
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new mrv_histogram object.
    #' @param ... If a param `Y` is given, the fit method is called with `...`.
    #' @return A new `mrv_histogram` object.
	initialize = function(...)
	{
		self$n_features = 0
		self$law_ = list()
		kwargs = list(...)
		if( !is.null(kwargs[["Y"]]) )
			base::do.call( self$fit , kwargs )
	},
	##}}}
	
	## fit ##{{{
	#' @description
    #' Fit method for the histograms
    #' @param Y [vector] Dataset to infer the histogram
    #' @param bins [list or vector or integer] bins values
    #' @return `self`
	fit = function( Y , bins = as.integer(100) )
	{
		if( !is.matrix(Y) ) Y = matrix( Y , nrow = length(Y) , ncol = 1 )
		self$n_features = ncol(Y)
		
		lbins = list()
		if( !is.list(bins) )
		{
			for( i in 1:self$n_features )
			{
				lbins[[i]] = bins
			}
		}
		else
		{
			lbins = bins
		}
		for( i in 1:self$n_features )
		{
			self$law_[[i]] = rv_histogram$new( Y = Y[,i] , lbins[[i]] )
		}
		return(self)
	},
	##}}}
	
	## rvs ##{{{
	#' @description
    #' Generation sample from the histogram
    #' @param n [integer] Number of samples drawn
    #' @return A matrix of samples
	rvs = function( n = 1 )
	{
		out = base::c()
		for( i in 1:self$n_features )
		{
			out =  base::rbind( out , self$law_[[i]]$rvs(n) )
		}
		return(base::t(out))
	},
	##}}}
	
	## cdf ##{{{
	#' @description
    #' Cumulative Distribution Function
    #' @param q [vector] Quantiles to compute the CDF
    #' @return cdf values
	cdf = function( q )
	{
		if( !is.matrix(q) ) q = matrix( q , nrow = length(q) , ncol = 1 )
		out = base::c()
		for( i in 1:self$n_features )
		{
			out =  base::rbind( out , self$law_[[i]]$cdf(q[,i]) )
		}
		return(base::t(out))
	},
	##}}}
	
	## sf ##{{{
	#' @description
    #' Survival Function
    #' @param q [vector] Quantiles to compute the SF
    #' @return sf values
	sf = function( q )
	{
		if( !is.matrix(q) ) q = matrix( q , nrow = length(q) , ncol = 1 )
		out = base::c()
		for( i in 1:self$n_features )
		{
			out =  base::rbind( out , self$law_[[i]]$sf(q[,i]) )
		}
		return(base::t(out))
	},
	##}}}
	
	## icdf ##{{{
	#' @description
    #' Inverse of Cumulative Distribution Function
    #' @param p [vector] Probabilities to compute the CDF
    #' @return icdf values
	icdf = function( p )
	{
		if( !is.matrix(p) ) p = matrix( p , nrow = length(p) , ncol = 1 )
		out = base::c()
		for( i in 1:self$n_features )
		{
			out =  base::rbind( out , self$law_[[i]]$icdf(p[,i]) )
		}
		return(base::t(out))
	},
	##}}}
	
	## isf ##{{{
	#' @description
    #' Inverse of Survival Function
    #' @param p [vector] Probabilities to compute the SF
    #' @return isf values
	isf = function( p )
	{
		if( !is.matrix(p) ) p = matrix( p , nrow = length(p) , ncol = 1 )
		out = base::c()
		for( i in 1:self$n_features )
		{
			out =  base::rbind( out , self$law_[[i]]$isf(p[,i]) )
		}
		return(base::t(out))
	}
	##}}}
	
	)
	
)
