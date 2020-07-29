
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
## Normal Law
##
################################################################################################

rv_GPD = R6::R6Class( "rv_GPD" ,
	
	inherit = RPOOstats::rv_abstract,
	
	############
	## Public ##
	############
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	loc   = 0,
	scale = 1,
	shape = 0,
	n_params = 3,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( loc = 0 , scale = 1 , shape = 0 )
	{
		super$initialize()
		self$loc   = loc
		self$scale = scale
		self$shape = shape
	},
	
	
	get_params = function()
	{
		invisible( base::c( self$loc , self$scale , self$shape ) )
	},
	
	
	set_params = function( params )
	{
		self$loc   = params[1]
		self$scale = params[2]
		self$shape = params[3]
	},
	
	
	#############
	## Methods ##
	#############
	
#	gradient_negloglikelihood = function( X )
#	{
#	},
	
	
	fit = function( X , method = "MLE" )
	{
		self$loc   = base::min(X)
		self$scale = 1
		self$shape = 0.1
		method = "Nelder-Mead"
#		method = "BFGS"
		stats::optim( par = base::c( self$loc , self$scale , self$shape ) , fn = private$optim_fct , X = X , method = method )
		
		invisible( list( loc = self$loc , scale = self$scale , shape = self$shape ) )
	},
	
	
	density = function( X )
	{
		Z = (X - self$loc) / self$scale
		if( self$shape == 0 )
		{
			return( base::exp( - Z ) / self$scale )
		}
		else
		{
			return( ( 1 + self$shape * Z )^( - 1 / self$shape - 1 ) / self$scale )
		}
	},
	
	
	rvs = function( size = 1 )
	{
		q = stats::runif( n = size , min = 0 , max = 1 )
		return( self$icdf(q) )
	},
	
	
	cdf = function( X )
	{
		Z = (X-self$loc) / self$scale
		Z[ Z < 0 ] = NaN
		if( self$shape == 0 )
		{
			return( 1 - base::exp(-Z) )
		}
		else
		{
			return( 1 - ( 1 + self$shape * Z )^( - 1 / self$shape ) )
		}
	},
	
	
	icdf = function( q )
	{
		if( self$shape == 0 )
		{
			return( self$loc - self$scale * base::log(1-q) )
		}
		else
		{
			return( self$loc + ( (1-q)^(-self$shape) - 1 ) * self$scale / self$shape )
		}
	},
	
	
	sf = function( X )
	{
		return( 1 - self$cdf(X) )
	},
	
	
	isf = function( q )
	{
		return( self$icdf(1-q) )
	}
	
	),
	
	private = list(
	
	#############
	## Methods ##
	#############
	
	
	)
)


