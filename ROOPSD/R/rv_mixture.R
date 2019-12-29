
##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2018                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the ARyga library. This  ##
## library makes it possible to study dynamic systems and to statistically      ##
## correct uni / multivariate data by applying optimal transport to             ##
## sparse histograms.                                                           ##
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
## that may rate  that it is complicated to manipulate,  and  that  also        ##
## therefore rates  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this rates that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2018                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## ARyga. Cette librairie permet d'étudier les systèmes dynamique et de         ##
## corriger statistiquement des données en uni/multivarié en appliquant le      ##
## transport optimal à des histogrammes creux.                                  ##
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


###############
## Libraries ##
###############

library(R6)


################################################################################################
##
## Exponential Law
##
################################################################################################

rv_mixture = R6::R6Class( "rv_mixture" ,
	
	inherit = RPOOstats::rv_abstract,
	
	############
	## Public ##
	############
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	n_laws = NULL,
	rv_laws = NULL,
	weights = NULL,
	low = NULL,
	high = NULL,
	n_params = 0,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( rv_laws , params = NULL , weights = base::rep( 1. / length(rv_laws) , length(rv_laws) ) )
	{
		super$initialize()
		self$n_laws = length(rv_laws)
		for( i in 1:self$n_laws )
		{
			self$rv_laws[[i]] = if( base::any( class(rv_laws[[i]]) %in% "rv_abstract" ) ) rv_laws[[i]] else base::do.call( rv_laws[[i]]$new , params[[i]] )
			self$n_params = self$n_params + self$rv_laws[[i]]$n_params
		}
		self$weights = weights / base::sum(weights)
		private$init_lhc()
	},
	
	
	###############
	## Accessors ##
	###############
	
	get_params = function()
	{
		params = base::c()
		for( i in 1:self$n_laws )
		{
			params = base::c( params , self$rv_laws[[i]]$get_params() )
		}
		
		params = base::c( params , self$weights )
		
		return(params)
	},
	
	
	set_params = function( params )
	{
		for( i in 1:self$n_laws )
		{
			self$rv_laws[[i]]$set_params( params[1:self$rv_laws[[i]]$n_params] )
			params = params[(self$rv_laws[[i]]$n_params+1):length(params)]
		}
		self$weights = params / base::sum(params)
	},
	
	
	#############
	## Methods ##
	#############
	
#	gradient_negloglikelihood = function( X )
#	{
#		dX = base::c()
#		
#		norm = self$density(X)
#		
#		## Gradient of each density
#		for( i in 1:self$n_laws )
#		{
#			dX = self$weights[i] * base::c( dX , self$rv_laws[[i]]$gradient_negloglikelihood(X) ) / norm
#		}
#		
#		
#		## Gradient of each weight
#		dX = base::c( 
#		
#		invisible(dX)
#	},
	
	
	fit = function( X , method = "MLE" , init = NULL , weights_init = NULL )
	{
		params_init = base::c()
		for( i in 1:self$n_laws )
		{
			if( is.null(init) )
			{
				self$rv_laws[[i]]$fit( X , method = "moments" )
			}
			else
			{
				self$rv_laws[[i]]$set_params( init[[i]] )
			}
			params_init = base::c( params_init , self$rv_laws[[i]]$get_params() )
		}
		
		
		weights_init = if( is.null(weights_init) ) base::rep( 1. / length(self$rv_laws) , length(self$rv_laws) ) else weights_init
		params_init = base::c( params_init , weights_init )
		
		stats::optim( params_init , private$optim_fct , X = X )
		
		private$init_lhc()
	},
	
	
	density = function( X )
	{
		density = matrix( 0 , ncol = length(X) , nrow = self$n_laws )
		for( i in 1:self$n_laws )
		{
			density[i,] = self$rv_laws[[i]]$density( X )
		}
		
		invisible( base::t(density) %*% self$weights )
	},
	
	
	rvs = function( size = 1 )
	{
		idx = base::sample( 1:self$n_laws , size = size , replace = TRUE , prob = self$weights )
		output = numeric(size)
		for( i in 1:self$n_laws )
		{
			idx_l = base::which( idx == i )
			output[idx_l] = self$rv_laws[[i]]$rvs( length(idx_l) )
		}
		
		invisible( output )
	},
	
	
	cdf = function( X )
	{
		cdf = matrix( 0 , ncol = length(X) , nrow = self$n_laws )
		for( i in 1:self$n_laws )
		{
			cdf[i,] = self$rv_laws[[i]]$cdf( X )
		}
		
		invisible( as.vector( base::t(cdf) %*% self$weights ) )
	},
	
	
	icdf = function( q )
	{
		invisible(private$icdffn(q))
	},
	
	
	sf = function( X )
	{
		invisible( 1. - self$cdf(X) )
	},
	
	
	isf = function( q )
	{
		invisible( self$icdf( 1. - q ) )
	}
	
	),
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	icdffn = NULL,
	
	
	#############
	## Methods ##
	#############
	
	init_lhc = function()
	{
		self$low  = .Machine$double.xmax
		self$high = -.Machine$double.xmax
		for( i in 1:self$n_laws )
		{
			self$low  = base::min( self$low  , self$rv_laws[[i]]$icdf(1e-6) )
			self$high = base::max( self$high , self$rv_laws[[i]]$icdf(1 - 1e-6) )
		}
		x = base::seq( self$low , self$high , length = 10000 )
		cdfx = self$cdf(x)
		private$icdffn = stats::approxfun( x = cdfx , y = x , yleft = 0 , yright = 1 )
	}
	
	)
)


