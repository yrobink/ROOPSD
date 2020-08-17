
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

base::rm( list = base::ls() )


###############
## Libraries ##
###############

library(R6)
library(devtools)
library(roxygen2)

roxygen2::roxygenize("../ROOPSD")
devtools::load_all("../ROOPSD")


###########################
## Useful plot functions ##
###########################

PlotTools = R6::R6Class( "PlotTools" , ##{{{
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	os = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
		self$os = self$get_os()
	},
	
	
	#############
	## Methods ##
	#############
	
	get_os = function()
	{
		sysinf = base::Sys.info()
		if( !is.null(sysinf) )
		{
			os = sysinf['sysname']
			if( os == 'Darwin' ) os = "osx"
		}
		else
		{
			## mystery machine
			os = .Platform$OS.type
			if( base::grepl( "^darwin"   , R.version$os ) ) os = "osx"
			if( base::grepl( "linux-gnu" , R.version$os ) ) os = "linux"
		}
		invisible(tolower(os))
	},
	
	new_screen = function( nrow , ncol , width = 5 * ncol , height = 5 * nrow )
	{
		if( self$os == "osx" )
		{
			grDevices::quartz( width = width , height = height )
		}
		if( self$os == "linux" )
		{
			grDevices::X11( width = width , height = height )
		}
		graphics::par( mfrow = base::c(nrow,ncol) )
	},
	
	wait = function()
	{
		while( base::names(grDevices::dev.cur()) !='null device' ) base::Sys.sleep(1)
	}
	
	)
)
##}}}

plt = PlotTools$new()


###############
## Functions ##
###############

## Tests
##======

test_parametric_law = function( claw , show = TRUE , ... ) ##{{{
{
	## Sample data
	params = list(...)
	law    = base::do.call( claw$new , params )
	X      = law$rvs(10000)
	
	## Fit
	lawf   = claw$new()
	lawf$fit(X)
	
	## cdf and sf
	x    = base::seq( base::min(X) , base::max(X) , length = 1000 )
	cdf  = law$cdf(x)
	sf   = law$sf(x)
	cdff = lawf$cdf(x)
	sff  = lawf$sf(x)
	
	## icdf and isf
	eps   = 1e-3
	q     = base::seq( eps , 1 - eps , length = 1000 )
	icdf  = law$icdf(q)
	isf   = law$isf(q)
	icdff = lawf$icdf(q)
	isff  = lawf$isf(q)
	
	## Density
	dens  = law$density(x)
	densf = lawf$density(x)
	hX    = hist( X , plot = FALSE )
	nb    = length(hX$breaks)
	hXp   = ( hX$breaks[2:nb] + hX$breaks[1:(nb-1)] ) / 2
	
	## Plot
	if( show )
	{
		plt$new_screen( nrow = 1 , ncol = 3 )
		
		## Plot cdf and sf
		plot(  x , cdf  , col = "red"  , type = "l" , main = law$name )
		lines( x , cdff , col = "blue" , type = "l" , lty = 2 )
		lines( x , sf   , col = "red"  , type = "l" )
		lines( x , sff  , col = "blue" , type = "l" , lty = 2 )
		
		## Plot icdf and isf
		plot(  x , icdf  , col = "red"  , type = "l" )
		lines( x , icdff , col = "blue" , type = "l" , lty = 2 )
		lines( x , isf   , col = "red"  , type = "l" )
		lines( x , isff  , col = "blue" , type = "l" , lty = 2 )
		
		## Plot density
		plot( hXp , hX$density , col = "red" , type = "h" , ylim = base::c(-0.01,base::max(dens,densf)) )
		lines( x , dens  , col = "red"  , type = "l" )
		lines( x , densf , col = "blue" , type = "l" , lty = 2 )
	}
	
}
##}}}

test_multivariate_normal = function( show = TRUE )##{{{
{
	lX    = list()
	xymin = 1e9
	xymax = -1e9
	for( i in 1:3 )
	{
		mean = stats::runif( n = 2 , min = -5 , max = 5 )
		cov  = ROOPSD::rspd_matrix( n = 1 , d = 2 , FALSE )
		lX[[i]] = ROOPSD::rmultivariate_normal( n = 10000 , mean = mean , cov = cov )
		xymin   = base::min( xymin , lX[[i]] )
		xymax   = base::max( xymax , lX[[i]] )
	}
	xylim = base::c(xymin,xymax)
	
	if( show )
	{
		plt$new_screen( nrow = 1 , ncol = 1 )
		
		color = base::c( "red" , "blue" , "green" )
		i = 1
		plot( lX[[i]][,1] , lX[[i]][,2] , xlim = xylim , ylim = xylim , col = color[i] , pch = 4 , xlab = "x" , ylab = "y" , main = "Bivariate Normal" )
		for( i in 2:3 )
			points( lX[[i]][,1] , lX[[i]][,2] , col = color[i] , pch = 4 )
		
	}
}
##}}}

test_rv_histogram = function( show = TRUE ) ##{{{
{
	## Sample data
	X = numeric(10000)
	X[1:2000] = stats::rexp( n = 2000 , rate = 1 )
	X[2001:10000] = stats::rnorm( n = 8000 , mean = 5 , sd = 1 )
	
	## Fit
	law = ROOPSD::rv_histogram$new( Y = X )
	
	## cdf and sf
	x    = base::seq( base::min(X) , base::max(X) , length = 1000 )
	cdf  = law$cdf(x)
	sf   = law$sf(x)
	
	## icdf and isf
	eps   = 1e-3
	q     = base::seq( eps , 1 - eps , length = 1000 )
	icdf  = law$icdf(q)
	isf   = law$isf(q)
	
	## Density
	dens  = law$density(x)
	hX    = hist( X , plot = FALSE )
	nb    = length(hX$breaks)
	hXp   = ( hX$breaks[2:nb] + hX$breaks[1:(nb-1)] ) / 2
	
	## Plot
	if( show )
	{
		plt$new_screen( nrow = 1 , ncol = 3 )
		
		## Plot cdf and sf
		plot(  x , cdf  , col = "red"  , type = "l" , main = law$name )
		lines( x , sf   , col = "red"  , type = "l" )
		
		## Plot icdf and isf
		plot(  x , icdf  , col = "red"  , type = "l" )
		lines( x , isf   , col = "red"  , type = "l" )
		
		## Plot density
		plot( hXp , hX$density , col = "red" , type = "h" , ylim = base::c(-0.01,base::max(dens)) )
		lines( x , dens  , col = "red"  , type = "l" )
	}
	
}
##}}}

test_rv_ratio_histogram = function( show = TRUE ) ##{{{
{
	## Sample data
	X = numeric(10000)
	X[1:2000] = 0
	X[2001:10000] = stats::rnorm( n = 8000 , mean = 5 , sd = 1 )
	
	## Fit
	law = ROOPSD::rv_ratio_histogram$new( Y = X , x0 = 0 )
	
	## cdf and sf
	x    = base::seq( base::min(X) , base::max(X) , length = 1000 )
	cdf  = law$cdf(x)
	sf   = law$sf(x)
	
	## icdf and isf
	eps   = 1e-3
	q     = base::seq( eps , 1 - eps , length = 1000 )
	icdf  = law$icdf(q)
	isf   = law$isf(q)
	
	## Plot
	if( show )
	{
		plt$new_screen( nrow = 1 , ncol = 2 )
		
		## Plot cdf and sf
		plot(  x , cdf  , col = "red"  , type = "l" , main = law$name )
		lines( x , sf   , col = "red"  , type = "l" )
		
		## Plot icdf and isf
		plot(  x , icdf  , col = "red"  , type = "l" )
		lines( x , isf   , col = "red"  , type = "l" )
		
	}
	
}
##}}}

## All in one
##===========

run_all_tests = function( show = FALSE )##{{{
{
	## Test parametric laws
	test_parametric_law( ROOPSD::Uniform     , show = show , min = -2 , max = 1 )
	test_parametric_law( ROOPSD::Normal      , show = show , mean = 3 , sd = 1.5 )
	test_parametric_law( ROOPSD::Exponential , show = show , rate = 0.5 )
	test_parametric_law( ROOPSD::Gamma       , show = show , scale = 1.5 , shape = 0.5 )
	test_parametric_law( ROOPSD::GEV         , show = show , loc = 1 , scale = 0.5 , shape = -0.6 )
	test_parametric_law( ROOPSD::GPD         , show = show , loc = 1 , scale = 0.5 , shape = -0.3 )
	test_multivariate_normal( show = show )
	
	## Non parametric
	test_rv_histogram( show = show )
	test_rv_ratio_histogram( show = show )
}
##}}}


##########
## main ##
##########


## Read command line arguments and run (or not) tests
##================================================{{{

args = commandArgs( trailingOnly = TRUE )
args_verbose = FALSE
args_run     = FALSE
if( length(args) > 0 )
{
	for( a in args )
	{
		if( a == "-r" || a == "--run-all-tests" )
			args_run = TRUE
		if( a == "-v" || a == "--verbose" )
			args_verbose = TRUE
	}
}

if( args_run )
	run_all_tests(args_verbose)

##}}}

plt$wait()
print("Done")

