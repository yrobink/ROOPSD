
################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2020, 2023                                          ##
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


#' GPD 
#'
#' @description
#' GPD distribution in OOP way. Based on AbstractDist
#'
#' @details
#' See AbstractDist for generic methods
#'
#' @examples
#' ## Generate sample
#' loc   = 0
#' scale = 0.5
#' shape = -0.3
#' gpd = ROOPSD::GPD$new( loc = loc , scale = scale , shape = shape )
#' X   = gpd$rvs( n = 1000 )
#'
#' ## And fit parameters
#' gpd$fit( X , loc = 0 )
#'
#' @export
GPD = R6::R6Class( "GPD",
	
	inherit = AbstractDist,
	
	## Private elements
	##==============={{{
	private = list(
	
	## Arguments
	##==========
	
	#' @field loc [double] location of the GPD law, fixed
	.loc = NULL,
	#' @field scale [double] scale of the GPD law
	.scale = NULL,
	#' @field shape [double] shape of the GPD law
	.shape = NULL,
	#' @field params [vector] params of the GPD law
	.params = NULL,
	
	## Methods
	##========
	
	fit_initialization = function(Y)##{{{
	{
		lmom = Lmoments::Lmoments(Y[Y>self$loc])
		
		itau  = lmom[1] / lmom[2]
		
		self$scale = lmom[1] * ( itau - 1 )
		self$shape = 2 - itau
	},
	##}}}
	
	gradient_negloglikelihood = function( params , Y )##{{{
	{
		self$params = params
		
		## Remove 0 from shape
		shape = self$shape
		shape[base::abs(shape) < 1e-10] = 1e-10
		
		## Usefull values
		Yu = Y[ Y > self$loc ]
		Z        = ( Yu - self$loc ) / self$scale
		ZZ       = 1. + shape * Z
		exponent = 1. + 1. / shape
		
		## Gradient
		dp    = base::c(NA,NA)
		dp[1] = base::sum( - exponent * shape * Z / ZZ / self$scale + 1. / self$scale )
		dp[2] = base::sum( - base::log(ZZ) / shape^2 + exponent * Z / ZZ )
		return(dp)
	}
	##}}}
	
	),
	##}}}
	
	## Active elements
	##================{{{
	active = list(
	
	## params ##{{{
	params = function(value)
	{
		if(missing(value))
		{
			return( list( loc = private$.loc , scale = private$.scale , shape = private$.shape ) )
		}
		else
		{
			if(is.numeric(value) && length(value) == 3 )
			{
				if( value[2] > 0 )
					private$.scale = value[2]
				private$.loc   = value[1]
				private$.shape = value[3]
			}
			
		}
	},
	##}}}
	
	## loc ##{{{
	loc = function(value)
	{
		if(missing(value))
		{
			return(private$.loc)
		}
		else
		{
			private$.loc = value
		}
	},
	##}}}
	
	## scale ##{{{
	scale = function(value)
	{
		if(missing(value))
		{
			return(private$.scale)
		}
		else
		{
			if(value > 0)
				private$.scale = value
		}
	},
	##}}}
	
	## shape ##{{{
	shape = function(value)
	{
		if(missing(value))
		{
			return(private$.shape)
		}
		else
		{
			private$.shape = value
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
	
	## Constructor
	##============
	
	## initialize ##{{{
	#' @description
    #' Create a new GPD object.
	#' @param loc   [double] location parameter
	#' @param scale [double] scale parameter
	#' @param shape [double] shape parameter
	#' @return A new `GPD` object.
	initialize = function( loc = 0 , scale = 1 , shape = -0.1 )
	{
		super$initialize( ROOPSD::dgpd , ROOPSD::pgpd , ROOPSD::qgpd , ROOPSD::rgpd , "GPD" , TRUE )
		self$loc   = loc
		self$scale = scale
		self$shape = shape
	},
	##}}}
	
	## Fit
	##====
	
	## fit ##{{{
	#' @description
    #' Fit method
    #' @param Y [vector] Dataset to infer the histogram
    #' @param loc [double] location parameter, if NULL used min(Y)
    #' @return `self`
	fit = function( Y , loc = NULL )
	{
		self$loc = 0
		if( is.null(loc) )
			loc = base::min(Y)
		super$fit( Y - loc )
		self$loc = loc
		
		return(self)
	}
	##}}}
	
	
	)
	##}}}
	
)



