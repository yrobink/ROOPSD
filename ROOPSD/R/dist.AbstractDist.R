
################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2020-2022                                           ##
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
## Copyright Yoann Robin, 2020-2022                                           ##
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
#' @description
#' Base class for OOP statistical distribution
#'
#' @details
#' This class is only used to be herited
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#' @importFrom numDeriv grad
#'
#' @export
AbstractDist = R6::R6Class( "AbstractDist",
	
	## Private elements
	##==============={{{
	
	private = list(
	
	## Arguments
	##==========
	
	#' @field name [string] name of the distribution
	.name = NULL,
	.has_gr_nlll = FALSE,
	#' @field opt [stats::optim result] Result of the MLE to find parameters
	.opt = NULL,
	#' @field cov [matrix] Covariance matrix of parameters, inverse of hessian
	.cov = NULL,
	#' @field coef [vector] Vector of coefficients
	.coef = NULL,
	
	## Methods
	##========
	
	## negloglikelihood ##{{{
	negloglikelihood = function( params , Y )
	{
		self$params = params
		return( - base::sum( base::suppressWarnings(self$logdensity(Y)) ) )
	},
	##}}}
	
	## gradient_negloglikelihood_withoutwarnings ##{{{
	gradient_neglll_ww = function( params , Y )
	{
		return( base::suppressWarnings(private$gradient_negloglikelihood( params , Y )) )
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
	name = function(value) 
	{
		if(base::missing(value))
		{
			return(private$.name)
		}
	},
	##}}}
	
	## name ##{{{
	coef = function(coef) 
	{
		if(base::missing(coef))
		{
			coef = c()
			for( i in 1:length(self$params) )
				coef = base::c(coef,self$params[[i]])
			return(coef)
		}
	},
	##}}}
	
	## opt ##{{{
	opt = function(value) 
	{
		if(base::missing(value))
		{
			return(private$.opt)
		}
	},
	##}}}
	
	## cov ##{{{
	cov = function(value) 
	{
		if(base::missing(value))
		{
			return(private$.cov)
		}
		else
		{
			private$.cov = value
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
	#' @field ddist [function] density function
	ddist = NULL,
	#' @field pdist [function] distribution function
	pdist = NULL,
	#' @field qdist [function] quantile function
	qdist = NULL,
	#' @field rdist [function] random generator function
	rdist = NULL,
	#' @field ks.test [ks.test] Goodness of fit with ks.test
	ks.test = NULL,
	#' @field fit_success [bool] TRUE only if the fit is a success and is occurred
	fit_success = FALSE,
	
	
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
	#' @param has_gr_nlll [bool] If the derived class has defined the gradient
	#'                           of the negative log-likelihood
	#' @return A new `AbstractDist` object.
	initialize = function( ddist , pdist , qdist , rdist , name , has_gr_nlll )
	{
		self$ddist = ddist
		self$pdist = pdist
		self$qdist = qdist
		self$rdist = rdist
		private$.name = name
		private$.has_gr_nlll = has_gr_nlll
	},
	##}}}
	
	## Methods
	##========
	
	## rvs ##{{{
	#' @description
    #' Generation sample from the histogram
    #' @param n [integer] Number of samples drawn
    #' @return [vector] A vector of samples
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
    #' @param x [vector] Values to compute the density
    #' @return [vector] density
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
    #' @param x [vector] Values to compute the log-density
    #' @return [vector] log of density
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
    #' @param q [vector] Quantiles to compute the CDF
    #' @return [vector] cdf values
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
    #' @param q [vector] Quantiles to compute the SF
    #' @return [vector] sf values
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
    #' @param p [vector] Probabilities to compute the CDF
    #' @return [vector] icdf values
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
    #' @param p [vector] Probabilities to compute the SF
    #' @return [vector] isf values
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
    #' @param Y [vector] Dataset to infer the histogram
    #' @param n_max_try [integer] Because the optim function can fails, the fit
    #'                            is retry n_try times.
    #' @return `self`
	fit = function( Y , n_max_try = 100 )
	{
		## Initialization
		private$fit_initialization(Y)
		
		## Prepare optimization params
		optparams = list( fn = private$negloglikelihood , method = "BFGS" , hessian = TRUE , Y = Y )
		if( private$.has_gr_nlll )
			optparams$gr = private$gradient_neglll_ww
		
		## Prepare for random initial condition
		params_m = self$coef
		params_c = diag(length(self$params)) / 10
		
		## Loop for the fit
		private$.opt = try( stop() , silent = TRUE )
		n_try = 0
		while( is(private$.opt,"try-error") && n_try < n_max_try )
		{
			## Try the fit
			n_try         = n_try + 1
			optparams$par = as.vector(self$params)
			private$.opt  = try( do.call( stats::optim , optparams ) , silent = TRUE )
			
			## Fail, draw a new initial condition
			if( is(private$.opt,"try-error") )
			{
				self$params = rmultivariate_normal( n = 1 , mean = params_m , cov = params_c )
				
				if( n_try %% 10 == 0 )
					params_c = 2 * params_c
			}
		}
		
		## OK, the fit is maybe really impossible
		if( is(private$.opt,"try-error") )
		{
			self$params = params_m
			return(self)
		}
		
		## Good!
		self$fit_success = TRUE
		self$params = self$opt$par
		self$cov    = base::try( base::solve(self$opt$hessian) , silent = TRUE )
		
		## Goodness of fit
		ksparams     = self$params
		ksparams$x   = Y
		ksparams$y   = self$pdist
		self$ks.test = suppressWarnings( base::do.call( stats::ks.test , ksparams ) )
		self$ks.test$data.name = NULL
		
		return(self)
	},
	##}}}
	
	## Confidence interval and diagnostic
	
	## qgradient ##{{{
	#' @description
	#' Gradient of the quantile function
    #' @param p [vector] Probabilities
    #' @param lower.tail [bool] If CDF or SF.
    #' @return [vector] gradient
	qgradient = function( p , lower.tail = TRUE )
	{
		if(!lower.tail)
			p = 1 - p
		
		FUN = function(c,p)
		{
			lc = as.list(c)
			names(lc) = names(self$params)
			lc$p = p
			return(do.call( self$qdist , lc ) )
		}
		g = matrix( NA , ncol = length(self$params) , nrow = length(p) )
		for( i in 1:length(p) )
		{
			g[i,] = numDeriv::grad( FUN , unlist(self$params) , p = p[i] )
		}
		
		return(g)
	},
	##}}}
	
	## qdeltaCI ##{{{
	#' @description
	#' Confidence interval of the quantile function
    #' @param p [vector] Probabilities
    #' @param Rt [bool] if Probabilities or return times
    #' @param alpha [double] level of confidence interval
    #' @return [list] Quantiles, and confidence interval
	qdeltaCI = function( p , Rt = FALSE , alpha = 0.05 )
	{
		if(Rt) p = 1. / p
		grad = as.matrix( self$qgradient( p , FALSE ) , nrow = length(p) )
		qvar = colSums(base::t(grad) * (self$cov %*% base::t(grad)))
		names(qvar) = base::c()
		
		qlevel       = self$isf(p)
		qlevel_left  = qlevel + stats::qnorm(     alpha / 2 ) * base::sqrt(qvar)
		qlevel_right = qlevel + stats::qnorm( 1 - alpha / 2 ) * base::sqrt(qvar)
		
		return( list(q = qlevel , left = qlevel_left , right = qlevel_right ) )
	},
	##}}}
	
	## pdeltaCI ##{{{
	#' @description
	#' Confidence interval of the CDF function
    #' @param x [vector] Quantiles
    #' @param Rt [bool] if Probabilities or return times
    #' @param alpha [double] level of confidence interval
    #' @return [list] CDF, and confidence interval
	pdeltaCI = function( x , Rt = FALSE , alpha = 0.05 )
	{
		grad = as.matrix( self$pgradient( x , FALSE ) , nrow = length(x) )
		pvar = colSums(base::t(grad) * (self$cov %*% base::t(grad)))
		names(pvar) = base::c()
		
		plevel       = self$sf(x)
		plevel_left  = plevel + stats::qnorm(     alpha / 2 ) * base::sqrt(pvar)
		plevel_right = plevel + stats::qnorm( 1 - alpha / 2 ) * base::sqrt(pvar)
		plevel_left  = base::pmax( 0 , plevel_left )
		plevel_left  = base::pmin( 1 , plevel_left )
		plevel_right = base::pmax( 0 , plevel_right )
		plevel_right = base::pmin( 1 , plevel_right )
		
		if(Rt)
		{
			plevel = 1. / plevel
			left   = 1. / plevel_right
			right  = 1. / plevel_left
			plevel_left  = left
			plevel_right = right
		}
		
		return( list(cdf = plevel , left = plevel_left , right = plevel_right ) )
	},
	##}}}
	
	## diagnostic ##{{{
	#' @description
	#' Diagnostic of the fitted law
    #' @param Y [vector] data to check
    #' @param alpha [double] level of confidence interval
    #' @return [NULL]
	diagnostic = function( Y , alpha = 0.05 )
	{
		graphics::par( mfrow = base::c(2,2) )
		rvY = rv_histogram$new(Y = Y)
		
		## Probability plot
		emp = rvY$cdf(Y)
		mod = self$cdf(Y)
		plot( mod , emp , xlab = "Model" , ylab = "Empirical" , xlim = base::c(0,1) , ylim = base::c(0,1) , main = "Probability" )
		lines( c(0,1) , c(0,1) , col = "blue" )
		
		## Quantile plot
		p   = seq( 0.01 , 0.99 , length = 100 )
		emp = rvY$icdf(p)
		mod = self$icdf(p)
		xylim = base::c(min(emp,mod),max(emp,mod))
		plot( mod , emp , xlim = xylim , ylim = xylim , xlab = "Model" , ylab = "Empirical" , main = "Quantile ")
		lines( xylim , xylim , col = "blue" )
		
		## Return level plot
		p   = self$cdf(Y)
		Rts = 10^base::seq( base::log10( base::min(1/p) )  , base::log10( base::max(1/p) ) , length = 100 )
		p   = 1. / Rts
		qlevel = self$qdeltaCI( p , FALSE , alpha = alpha )
		ylim = base::c(base::min(qlevel$left,Y),base::max(qlevel$right,Y))
		plot(  Rts , qlevel$q     , col = "black" , type = "l" , log = "x" , main = "Return level" , xlab = "Return time" , ylab = "Level" , ylim = ylim )
		lines( Rts , qlevel$left  , col = "blue" )
		lines( Rts , qlevel$right , col = "blue" )
		points( 1. / rvY$sf(Y) , Y )
		
		## Density plot
		minY = base::min(Y)
		maxY = base::max(Y)
		delY = 0.1 * (maxY-minY)
		x = base::seq( minY - delY , maxY + delY , length = 1000 )
		hist( Y , breaks = min(floor(0.1*length(Y)) + 1,100) , col = grDevices::rgb(1,0,0,0.1) , freq = FALSE , xlab = "Value" , main = "Density" )
		lines( x , self$density(x) , col = "red" , type = "l" )
		points( Y , base::rep(0,length(Y)) )
		
		invisible(NULL)
	}
	##}}}
	
	)
	##}}}
	
)

