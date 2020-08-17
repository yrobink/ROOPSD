
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


#' rv_mixture 
#'
#' @description
#' rv_mixture distribution in OOP way. 
#'
#' @details
#' No fit allowed.
#'
#' @examples
#' ## Define the mixture
#' l_dist  = list( Exponential$new() , Normal$new( mean = 5 , sd = 1 ) )
#' weights = base::c( 0.2 , 0.8 )
#' rvX = rv_mixture$new( l_dist , weights )
#'
#' ## Draw samples
#' X   = rvX$rvs( 1000 )
#'
#' @export
rv_mixture = R6::R6Class( "rv_mixture" ,
	
	#####################
	## Active elements ##
	#####################
	
	active = list(
	
	## weights ##{{{
	weights = function(value) 
	{
		if(base::missing(value))
		{
			return(private$.weights)
		}
		else
		{
			if( is.null(value) )
				value = numeric(self$n_dist) + 1
			private$.weights = value / base::sum(value)
		}
	},
	##}}}
	
	## n_dist ##{{{
	n_dist = function(value) 
	{
		if(base::missing(value))
		{
			return(private$.n_dist)
		}
	},
	##}}}
	
	## l_dist ##{{{
	l_dist = function(value) 
	{
		if(base::missing(value))
		{
			return(private$.l_dist)
		}
	}
	##}}}
	
	),
	
	
	#####################
	## Public elements ##
	#####################
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field l_dist  [list]    List of distributions.
	#' @field n_dist  [integer] Numbers of distribution.
	#' @field weights [vector]  Weights of the distributions.
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new rv_mixture object.
	#' @param l_dist  [list]    List of ROOPSD distributions.
	#' @param weights [vector]  Weights of the distributions. If NULL,
	#'                          1 / length(l_dist) is used.
    #' @return A new `rv_mixture` object.
	initialize = function( l_dist , weights = NULL )
	{
		private$.l_dist = l_dist
		private$.n_dist = length(l_dist)
		self$weights    = weights
		
		## icdf function
		X = self$rvs(10000)
		xmin = base::min(X)
		xmax = base::max(X)
		dx   = 0.1 * (xmax - xmin)
		x    = base::seq( xmin , xmax , length = 1000 )
		cdf  = self$cdf(x)
		private$.icdf = stats::approxfun( cdf , x , yleft = 0 , yright = 1 , ties = "ordered" )
	},
	##}}}
	
	## rvs ##{{{
	#' @description
    #' Generation sample from the histogram
    #' @param n [integer] Number of samples drawn
    #' @return A vector of samples
	rvs = function( n )
	{
		## Build tools to apply rvs
		##=========================
		n_per_dist = as.integer(n * self$weights)
		n_per_dist[1] = n_per_dist[1] + n - base::sum(n_per_dist)
		
		rvsdist = list()
		for( i in 1:self$n_dist )
			rvsdist[[i]] = list( dist = self$l_dist[[i]] , n = n_per_dist[i] )
		
		## And draw values
		##================
		res = unlist( lapply( rvsdist , function(rvs) { rvs$dist$rvs(rvs$n) } ) )
		
		return(res)
	},
	##}}}
	
	## density ##{{{
	#' @description
    #' Density function
    #' @param x [vector] Values to compute the density
    #' @return density
	density = function( x ) 
	{
		y = 0 * x
		for( i in 1:self$n_dist )
		{
			y = y + self$l_dist[[i]]$density(x) * self$weights[i]
		}
		return(y)
	},
	##}}}
	
	## logdensity ##{{{
	#' @description
    #' Log density function
    #' @param x [vector] Values to compute the log-density
    #' @return the log density
	logdensity = function( x ) 
	{
		return(base::log(self$density(x)))
	},
	##}}}
	
	## cdf ##{{{
	#' @description
    #' Cumulative Distribution Function
    #' @param q [vector] Quantiles to compute the CDF
    #' @return cdf values
	cdf = function( q )
	{
		p = 0 * q
		for( i in 1:self$n_dist )
		{
			p = p + self$l_dist[[i]]$cdf(q) * self$weights[i]
		}
		return(p)
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
	
	.l_dist  = NULL,
	.n_dist  = NULL,
	.weights = NULL,
	.icdf    = NULL
	)
)


