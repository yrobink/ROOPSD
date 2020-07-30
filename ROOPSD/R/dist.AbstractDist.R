
################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2020                                                ##
##                                                                            ##
## yoann.robin.k@gmail.com                                                    ##
##                                                                            ##
## This software is a computer program that is part of the ROOPSD (R Object   ##
## Oriented Programming for Statistical Distribution). This library makes it  ##
## possible to perform bias to  draw, fit and access to CDF and SF function   ##
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

#' AbstractDist 
#'
#' Base class for OOP statistical distribution
#'
#' @export
AbstractDist = R6::R6Class( "AbstractDist",
	
	## Private elements
	##==============={{{
	
	private = list(
	
	## Arguments
	##==========
	
	#' @field name name of the distribution
	.name = NULL,
	
	## Methods
	##========
	
	## negloglikelihood ##{{{
	negloglikelihood = function( params , Y )
	{
		self$params = params
		return(-base::sum(self$logdensity(Y)))
	},
	##}}}
	
	## fit_initialization ##{{{
	fit_initialization = function(Y)
	{}
	##}}}
	
	),
	##}}}
	
	## Active elements
	##================{{{
	active = list(
	
	## name ##{{{
	#' @description
    #' Setter/getter of name, set is protected
	name = function(value) 
	{
		if(base::missing(value))
		{
			return(private$.name)
		}
	}
	##}}}
	
	),
	##}}}
	
	## Public elements
	##============={{{
	public = list(
	
	## Arguments
	##==========
	#' @field ddist density function
	ddist = NULL,
	#' @field pdist distribution function
	pdist = NULL,
	#' @field qdist quantile function
	qdist = NULL,
	#' @field rdist random generator function
	rdist = NULL,
	
	## Constructor
	##============
	
	## initialize ##{{{
	#' @description
    #' Create a new AbstractDist object.
	#' @param ddist [function] Density function, e.g. dnorm
	#' @param pdist [function] Distribution function, e.g. pnorm
	#' @param qdist [function] Quantile function, e.g. qnorm
	#' @param rdist [function] Random generator function, e.g. rnorm
	#' @param name  [str]      name of the distribution
	#' @return A new `AbstractDist` object.
	initialize = function( ddist , pdist , qdist , rdist , name )
	{
		self$ddist = ddist
		self$pdist = pdist
		self$qdist = qdist
		self$rdist = rdist
		private$.name = name
	},
	##}}}
	
	## Methods
	##========
	
	## rvs ##{{{
	#' @description
    #' Generation sample from the histogram
    #' @param n Number of samples drawn
    #' @return A vector of samples
	rvs = function( n )
	{
		params = self$params
		params$n = n
		return(base::do.call( self$rdist , params ))
	},
	##}}}
	
	## density ##{{{
	#' @description
    #' Density function
    #' @param x Values to compute the density
    #' @return density
	density = function(x)
	{
		params = self$params
		params$x = x
		return(base::do.call( self$ddist , params ) )
	},
	##}}}
	
	## logdensity ##{{{
	#' @description
    #' Log density function
    #' @param x Values to compute the log-density
    #' @return log of density
	logdensity = function(x)
	{
		params = self$params
		params$x = x
		params$log = TRUE
		return(base::do.call( self$ddist , params ) )
	},
	##}}}
	
	## cdf ##{{{
	#' @description
    #' Cumulative Distribution Function
    #' @param q Quantiles to compute the CDF
    #' @return cdf values
	cdf = function(q)
	{
		params = self$params
		params$q = q
		return(base::do.call( self$pdist , params ) )
	},
	##}}}
	
	## sf ##{{{
	#' @description
    #' Survival Function
    #' @param q Quantiles to compute the SF
    #' @return sf values
	sf = function(q)
	{
		params = self$params
		params$q = q
		params$lower.tail = FALSE
		return(base::do.call( self$pdist , params ) )
	},
	##}}}
	
	## icdf ##{{{
	#' @description
    #' Inverse of Cumulative Distribution Function
    #' @param p Probabilities to compute the CDF
    #' @return icdf values
	icdf = function(p)
	{
		params = self$params
		params$p = p
		params$lower.tail = TRUE
		return(base::do.call( self$qdist , params ) )
	},
	##}}}
	
	## isf ##{{{
	#' @description
    #' Inverse of Survival Function
    #' @param p Probabilities to compute the SF
    #' @return isf values
	isf = function(p)
	{
		params = self$params
		params$p = p
		params$lower.tail = FALSE
		return(base::do.call( self$qdist , params ) )
	},
	##}}}
	
	## Fit
	##====
	
	## fit ##{{{
	#' @description
    #' Fit method
    #' @param Y Dataset to infer the histogram
    #' @return `self`
	fit = function(Y)
	{
		private$fit_initialization(Y)
		opt = stats::optim( par = as.vector(self$params) , fn = private$negloglikelihood , gr = private$gradient_negloglikelihood , method = "BFGS" , Y = Y )
		self$params = opt$par
		
		return(self)
	}
	##}}}
	
	)
	##}}}
	
)

