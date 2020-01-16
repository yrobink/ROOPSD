
##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2020                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the ROOPSD (R Object     ##
## Oriented Programming for Statistical Distribution). This library makes it    ##
## possible to perform bias to  draw, fit and access to CDF and SF function for ##
## for classic statistical distribution.                                        ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2020                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## ROOPSD (R Object Oriented Programming for Statistical Distribution). Cette   ##
## librarie permet de tirer aléatoirement, inférer et acceder aux fonctions CDF ##
## et SF des distributions statistisques classiques.                            ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################

#' AbstractDist 
#'
#' Base class for OOP statistical distribution
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param ddist [function]
#'        Density function, e.g. dnorm
#' @param pdist [function]
#'        Distribution function, e.g. pnorm
#' @param qdist [function]
#'        Quantile function, e.g. qnorm
#' @param rdist [function]
#'        Random generator function, e.g. rnorm
#' @param freeze [boolean]
#'        If we freeze or note the distribution. TRUE after fit in general
#' @param n  [integer]
#'        For rvs function, number of samples drawn
#' @param x  [integer]
#'        For density and logdensity function. Vector of quantiles
#' @param q  [integer]
#'        For cdf and sf functions. Vector of quantiles
#' @param p  [integer]
#'        For icdf and isf functions. Vector of probabilities
#' @param Y  [integer]
#'        For fit function. Dataset to infer parameters. (max likelihood)
#'
#' @return Object of \code{\link{R6Class}} with methods to use a statistical distribution.
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(ddist,pdist,qdist,rdist,freeze)}}{This method is used to create object of this class with \code{AbstractDist}}
#'   \item{\code{rvs(n)}}{Draw n samples}.
#'   \item{\code{density(x)}}{Density along vector of quantile x}.
#'   \item{\code{logdensity(x)}}{log of density along vector of quantile x}.
#'   \item{\code{cdf(q)}}{Cumulative Distribution Function along vector of quantile q}.
#'   \item{\code{sf(q)}}{Survival function (1-CDF) along vector of quantile q}.
#'   \item{\code{icdf(p)}}{Inverse of cdf along vector of probabilities p}.
#'   \item{\code{isf(p)}}{Inverse of sf along vector of probabilities p}.
#'   \item{\code{fit(Y)}}{Fit function to infer parameters. By max likelihood}.
#' }
#' @examples
#' ##
#' ##
#' @export
AbstractDist = R6::R6Class( "AbstractDist",
	
	## Public elements
	##============={{{
	public = list(
	
	## Arguments
	##==========
	ddist = NULL,
	pdist = NULL,
	qdist = NULL,
	rdist = NULL,
	freeze = FALSE,
	
	## Constructor
	##============
	
	initialize = function( ddist , pdist , qdist , rdist , freeze = FALSE )##{{{
	{
		self$ddist = ddist
		self$pdist = pdist
		self$qdist = qdist
		self$rdist = rdist
		self$freeze = freeze
	},
	##}}}
	
	## Methods
	##========
	
	rvs = function( n )##{{{
	{
		params = private$params()
		params$n = n
		return(base::do.call( self$rdist , params ))
	},
	##}}}
	
	density = function(x)##{{{
	{
		params = private$params()
		params$x = x
		return(base::do.call( self$ddist , params ) )
	},
	##}}}
	
	logdensity = function(x)##{{{
	{
		params = private$params()
		params$x = x
		params$log = TRUE
		return(base::do.call( self$ddist , params ) )
	},
	##}}}
	
	cdf = function(q)##{{{
	{
		params = private$params()
		params$q = q
		return(base::do.call( self$pdist , params ) )
	},
	##}}}
	
	sf = function(q)##{{{
	{
		params = private$params()
		params$q = q
		params$lower.tail = FALSE
		return(base::do.call( self$pdist , params ) )
	},
	##}}}
	
	icdf = function(p)##{{{
	{
		params = private$params()
		params$p = p
		params$lower.tail = TRUE
		return(base::do.call( self$qdist , params ) )
	},
	##}}}
	
	isf = function(p)##{{{
	{
		params = private$params()
		params$p = p
		params$lower.tail = FALSE
		return(base::do.call( self$qdist , params ) )
	},
	##}}}
	
	## Fit
	##====
	
	fit = function(Y)##{{{
	{
		private$fit_initialization(Y)
		opt = stats::optim( par = as.vector(private$params()) , fn = private$negloglikelihood , method = "BFGS" , Y = Y )
		private$set_params(opt$par)
	}
	##}}}
	
	),
	##}}}
	
	## Private elements
	##==============={{{
	
	private = list(
	
	params = function()##{{{
	{},
	##}}}
	
	set_params = function(params)##{{{
	{},
	##}}}
	
	negloglikelihood = function( params , Y )##{{{
	{
		private$set_params(params)
		return(-base::sum(self$logdensity(Y)))
	},
	##}}}
	
	fit_initialization = function(Y)##{{{
	{}
	##}}}
	
	)
	##}}}
	
)

