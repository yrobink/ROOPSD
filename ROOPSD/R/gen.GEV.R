
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


## dgev {{{

#' dgev
#'
#' Density function of Generalized Extreme Value distribution
#'
#' @param x      [vector] Vector of values
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param log    [bool]   Return log of density if TRUE, default is FALSE
#'
#' @return [vector] Density of GEV at x
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' x = base::seq( -5 , 5 , length = 1000 )
#' y = dgev( x , loc = loc , scale = scale , shape = shape )
#' @export
dgev = function( x , loc = 0 , scale = 1 , shape = 0 , log = FALSE )
{
	size_x = length(x)
	loc   = if( length(loc)   == size_x ) loc   else base::rep( loc[1]   , size_x ) 
	scale = if( length(scale) == size_x ) scale else base::rep( scale[1] , size_x ) 
	shape = if( length(shape) == size_x ) shape else base::rep( shape[1] , size_x ) 
	
	
	Z     = ( x - loc ) / scale
	valid = (1 + shape * Z > 0)
	shape_zero  = ( base::abs(shape) < 1e-10 )
	cshape_zero = !shape_zero
	
	TX = numeric(length(x)) + NA
	
	if( base::any(shape_zero) )
	{
		TX[shape_zero] = base::exp( - Z[shape_zero] )
	}
	if( base::any(cshape_zero) )
	{
		TX[cshape_zero] = ( 1 + shape[cshape_zero] * Z[cshape_zero] )^( - 1. / shape[cshape_zero] )
	}
	
	out = TX^( shape + 1  ) * base::exp( - TX ) / scale
	if( base::any(!valid) )
	{
		out[!valid] = 0
	}
	
	if( log )
		return(base::log(out))
	else
		return(out)

}
##}}}

## pgev {{{

#' pgev
#'
#' Cumulative distribution function (or survival function) of Generalized 
#' Extreme Value distribution
#'
#' @param q      [vector] Vector of quantiles
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param lower.tail [bool] Return CDF if TRUE, else return survival function
#'
#' @return [vector] CDF (or SF) of GEV at x
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' x = base::seq( -5 , 5 , length = 1000 )
#' cdfx = pgev( x , loc = loc , scale = scale , shape = shape )
#' @export
pgev = function( q , loc = 0 , scale = 1 , shape = 0 , lower.tail = TRUE )
{
	if( !lower.tail )
	{
		return( 1. - pgev( q , loc , scale , shape , lower.tail = TRUE ) )
	}
	
	size_q = base::length(q)
	loc   = if( length(loc)   == size_q ) loc   else base::rep( loc[1]   , size_q )
	scale = if( length(scale) == size_q ) scale else base::rep( scale[1] , size_q )
	shape = if( length(shape) == size_q ) shape else base::rep( shape[1] , size_q )
	
	shape_zero  = ( base::abs(shape) < 1e-10 )
	cshape_zero = !shape_zero
	
	Z = ( q - loc ) / scale
	out = numeric(size_q) + NA
	
	if( base::any(shape_zero) )
	{
		out[shape_zero] = base::exp( - base::exp( - Z[shape_zero] ) )
	}
	
	if( base::any(cshape_zero) )
	{
		out[cshape_zero] = base::exp( - ( 1. + shape[cshape_zero] * Z[cshape_zero] )^( - 1. / shape[cshape_zero] ) )
	}
	
	valid = (1 + shape * Z > 0)
	if( base::any(!valid) )
	{
		out[(shape > 0) & !valid] = 0
		out[(shape < 0) & !valid] = 1
	}
	return(out)
}
##}}}

## qgev {{{

#' qgev
#'
#' Inverse of CDF (or SF) function of Generalized Extreme Value distribution
#'
#' @param p      [vector] Vector of probabilities
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param lower.tail [bool] Return inverse of CDF if TRUE, else return inverse 
#'                          of survival function
#'
#' @return [vector] Inverse of CDF or SF of GEV for probabilities p
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' p = base::seq( 0.01 , 0.99 , length = 100 )
#' q = qgev( p , loc = loc , scale = scale , shape = shape )
#' @export
qgev = function( p , loc = 0 , scale = 1 , shape = 0 , lower.tail = TRUE )
{
	if( !lower.tail )
	{
		return( qgev( 1. - p , loc , scale , shape , lower.tail = TRUE ) )
	}
	
	size_p = base::length(p)
	loc   = if( length(loc)   == size_p ) loc   else base::rep( loc[1]   , length(p) )
	scale = if( length(scale) == size_p ) scale else base::rep( scale[1] , length(p) )
	shape = if( length(shape) == size_p ) shape else base::rep( shape[1] , length(p) )
	
	shape_zero  = ( base::abs(shape) < 1e-10 )
	cshape_zero = !shape_zero
	
	out = numeric(length(p)) + NA
	if( base::any(shape_zero) )
	{
		out[shape_zero] = loc[shape_zero] - scale[shape_zero] * base::log( - base::log(p[shape_zero]) )
	}
	
	if( base::any(cshape_zero) )
	{
		out[cshape_zero] = loc[cshape_zero] + scale[cshape_zero] * ( ( - base::log(p[cshape_zero]) )^(- shape[cshape_zero]) - 1. ) / shape[cshape_zero]
	}
	
	return(out)
}
##}}}

## rgev {{{

#' rgev
#'
#' Random value generator of Generalized Extreme Value distribution
#'
#' @param n      [int]    Numbers of values generated
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#'
#' @return [vector] Random value following a GEV(loc,scale,shape)
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' gev = rgev( 100 , loc = loc , scale = scale , shape = shape )
#' @export
rgev = function( n = 1 , loc = 0 , scale = 1 , shape = 0 )
{
	p = stats::runif( n = n )
	return( qgev( p , loc , scale , shape ) )
}
##}}}



