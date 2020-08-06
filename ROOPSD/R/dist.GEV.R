
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


#' GEV 
#'
#' @description
#' GEV distribution in OOP way. Based on AbstractDist
#'
#' @details
#' See AbstractDist for generic methods
#'
#' @importFrom Lmoments Lmoments
#'
#' @examples
#' ## Generate sample
#' loc   = 0
#' scale = 0.5
#' shape = -0.3
#' gev = ROOPSD::GEV$new( loc = loc , scale = scale , shape = shape )
#' X   = gev$rvs( n = 1000 )
#'
#' ## And fit parameters
#' gev$fit(X)
#'
#' @export
GEV = R6::R6Class( "GEV",
	
	inherit = AbstractDist,
	
	## Private elements
	##==============={{{
	private = list(
	
	## Arguments
	##==========
	
	#' @field loc [double] location of the GEV law
	.loc = NULL,
	#' @field scale [double] scale of the GEV law
	.scale = NULL,
	#' @field shape [double] shape of the GEV law
	.shape = NULL,
	#' @field params [vector] params of the GEV law
	.params = NULL,
	
	## Methods
	##========
	
	fit_initialization = function(Y)##{{{
	{
		lmom = Lmoments::Lmoments(Y)
		
		tau3  = lmom[3] / lmom[2]
		co    = 2. / ( 3. + tau3 ) - base::log(2) / base::log(3)
		kappa = 7.8590 * co + 2.9554 * co^2
		g     = base::gamma( 1. + kappa )
		
		
		self$scale = lmom[2] * kappa / ( (1 - 2^( - kappa )) * g )
		self$loc   = lmom[1] - self$scale * (1 - g) / kappa
		self$shape = - kappa
	},
	##}}}
	
	gradient_negloglikelihood = function( params , Y )##{{{
	{
		self$params = params
		
		## Remove 0 from shape
		shape = self$shape
		shape[base::abs(shape) < 1e-10] = 1e-10
		
		## Usefull values
		Z      = ( Y - self$loc ) / self$scale
		Za1    = 1 + shape * Z
		ishape = 1. / shape
		Zamsi  = Za1^(-ishape)
		dp     = base::c(NA,NA,NA)
		
		
		## Test
		if( !(self$scale > 0) || !base::all(Za1 > 0) )
			return(dp)
		
		## Gradient
		dp[1] = base::sum( ( Zamsi - 1 - shape ) / ( self$scale * Za1 ) )
		dp[2] = base::sum( ( 1. + Z * ( Zamsi - 1 - shape ) / Za1 ) / self$scale )
		dp[3] = base::sum( ( ( Zamsi - 1. ) * base::log(Za1) * ishape**2 + ( 1. + ishape - ishape * Zamsi ) * Z / Za1 ) )
		
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
				private$.loc = value[1]
				if( value[2] > 0 )
					private$.scale = value[2]
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
    #' Create a new GEV object.
	#' @param loc   [double] location parameter
	#' @param scale [double] scale parameter
	#' @param shape [double] shape parameter
	#' @return A new `GEV` object.
	initialize = function( loc = 0 , scale = 1 , shape = -0.1 )
	{
		super$initialize( ROOPSD::dgev , ROOPSD::pgev , ROOPSD::qgev , ROOPSD::rgev , "GEV" , TRUE )
		self$loc   = loc
		self$scale = scale
		self$shape = shape
	}
	##}}}
	
	)
	##}}}
	
)


