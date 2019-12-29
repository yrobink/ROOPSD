
#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

rm( list = ls() )

###############
## Libraries ##
###############

library(devtools)
devtools::load_all("RPOOstats")
library(ggplot2)
library(gridExtra)
library(extRemes)


###############
## Functions ##
###############

test_uniform_law = function() ##{{{
{
	## Define law
	rvX = RPOOstats::rv_uniform$new( min = -2 , max = 3 )
	
	## Random number
	X = rvX$rvs(10000)
	
	## CDF
	x = base::seq( -2 , 3 , length = 1000 )
	cdfX = rvX$cdf(x)
	
	## iCDF
	q = base::seq( 0 , 1 , length = 1000 )
	icdfX = rvX$icdf(q)
	
	## Fit
	rvU = RPOOstats::rv_uniform$new()
	rvU$fit( X , method = "MLE" )
	print( rvU$get_params() )
	
	## Plot
	p1 = ggplot2::ggplot( data.frame( X = X ) , ggplot2::aes( x = X ) ) + ggplot2::geom_histogram( ggplot2::aes( y = ..density.. ) , binwidth = 0.05 , color = "red" , fill = "red" , alpha = 0.5 ) + ggplot2::geom_density( alpha = 0.2 , color = "red" , fill = "red" ) + ggplot2::ggtitle( "Histogram" )
	p2 = ggplot2::ggplot( data.frame( x = x , y = cdfX  ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "CDF" )
	p3 = ggplot2::ggplot( data.frame( x = q , y = icdfX ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "iCDF" )
	
	gridExtra::grid.arrange( p1 , p2 , p3 , nrow = 2 )

} ##}}}

test_normal_law = function() ##{{{
{
	## Define law
	rvX = RPOOstats::rv_normal$new( loc = 4 , scale = 3 )
	
	## Random number
	X = rvX$rvs(10000)
	
	## CDF
	x = base::seq( -6 , 14 , length = 1000 )
	cdfX = rvX$cdf(x)
	
	## iCDF
	q = base::seq( 0 , 1 , length = 1000 )
	icdfX = rvX$icdf(q)
	
	## Fit
	rvN = RPOOstats::rv_normal$new()
	rvN$fit( X , method = "MLE" )
	print( rvN$get_params() )
	
	## Plot
	p1 = ggplot2::ggplot( data.frame( X = X ) , ggplot2::aes( x = X ) ) + ggplot2::geom_histogram( ggplot2::aes( y = ..density.. ) , binwidth = 0.05 , color = "red" , fill = "red" , alpha = 0.5 ) + ggplot2::geom_density( alpha = 0.2 , color = "red" , fill = "red" ) + ggplot2::ggtitle( "Histogram" )
	p2 = ggplot2::ggplot( data.frame( x = x , y = cdfX  ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "CDF" )
	p3 = ggplot2::ggplot( data.frame( x = q , y = icdfX ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "iCDF" )
	
	gridExtra::grid.arrange( p1 , p2 , p3 , nrow = 2 )

} ##}}}

test_exponential_law = function() ##{{{
{
	## Define law
	rvX = RPOOstats::rv_exponential$new( rate = 3 )
	
	## Random number
	X = rvX$rvs(10000)
	
	## CDF
	x = base::seq( 0 , 5 , length = 1000 )
	cdfX = rvX$cdf(x)
	
	## iCDF
	q = base::seq( 0 , 1 , length = 1000 )
	icdfX = rvX$icdf(q)
	
	## Fit
	rvE = RPOOstats::rv_exponential$new()
	rvE$fit( X , method = "MLE" )
	print( rvE$get_params() )
	
	## Plot
	p1 = ggplot2::ggplot( data.frame( X = X ) , ggplot2::aes( x = X ) ) + ggplot2::geom_histogram( ggplot2::aes( y = ..density.. ) , binwidth = 0.05 , color = "red" , fill = "red" , alpha = 0.5 ) + ggplot2::geom_density( alpha = 0.2 , color = "red" , fill = "red" ) + ggplot2::ggtitle( "Histogram" )
	p2 = ggplot2::ggplot( data.frame( x = x , y = cdfX  ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "CDF" )
	p3 = ggplot2::ggplot( data.frame( x = q , y = icdfX ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "iCDF" )
	
	gridExtra::grid.arrange( p1 , p2 , p3 , nrow = 2 )

} ##}}}

test_gamma_law = function() ##{{{
{
	## Define law
	rvX = RPOOstats::rv_gamma$new( shape = 2 , rate = 3 )
	
	## Random number
	X = rvX$rvs(10000)
	
	## CDF
	x = base::seq( 0 , 5 , length = 1000 )
	cdfX = rvX$cdf(x)
	
	## iCDF
	q = base::seq( 0 , 1 , length = 1000 )
	icdfX = rvX$icdf(q)
	
	## Fit
	rvG = RPOOstats::rv_gamma$new()
	rvG$fit( X , method = "MLE" )
	print( rvG$get_params() )
	
	## Plot
	p1 = ggplot2::ggplot( data.frame( X = X ) , ggplot2::aes( x = X ) ) + ggplot2::geom_histogram( ggplot2::aes( y = ..density.. ) , binwidth = 0.05 , color = "red" , fill = "red" , alpha = 0.5 ) + ggplot2::geom_density( alpha = 0.2 , color = "red" , fill = "red" ) + ggplot2::ggtitle( "Histogram" )
	p2 = ggplot2::ggplot( data.frame( x = x , y = cdfX  ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "CDF" )
	p3 = ggplot2::ggplot( data.frame( x = q , y = icdfX ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "iCDF" )
	
	gridExtra::grid.arrange( p1 , p2 , p3 , nrow = 2 )

} ##}}}

test_mixture_law = function() ##{{{
{
	## Define law
	rvX = rv_mixture$new( list( rv_exponential$new(2) , rv_normal$new(3,0.5) , rv_normal$new(-3,0.2) ) , weights = base::c(0.2,0.7,0.1) )
	
	## Random number
	X = rvX$rvs(10000)
	
	## CDF
	x = base::seq( -5 , 5 , length = 1000 )
	cdfX = rvX$cdf(x)
	
	## iCDF
	q = base::seq( 0 , 1 , length = 1000 )
	icdfX = rvX$icdf(q)
	
	## Fit
	rvM = RPOOstats::rv_mixture$new( list( rv_exponential$new() , rv_normal$new() , rv_normal$new() ) )
	rvM$fit( X , method = "MLE" , init = list( c(1.5) , c(3.4,0.8) , c(-2,0.4) ) )
	print( rvM$get_params() )
	
	## Plot
	p1 = ggplot2::ggplot( data.frame( X = X ) , ggplot2::aes( x = X ) ) + ggplot2::geom_histogram( ggplot2::aes( y = ..density.. ) , binwidth = 0.05 , color = "red" , fill = "red" , alpha = 0.5 ) + ggplot2::geom_density( alpha = 0.2 , color = "red" , fill = "red" ) + ggplot2::ggtitle( "Histogram" )
	p2 = ggplot2::ggplot( data.frame( x = x , y = cdfX  ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "CDF" )
	p3 = ggplot2::ggplot( data.frame( x = q , y = icdfX ) , ggplot2::aes(x = x , y = y ) ) + ggplot2::geom_line( color = "red" ) + ggplot2::ggtitle( "iCDF" )
	
	gridExtra::grid.arrange( p1 , p2 , p3 , nrow = 2 )

} ##}}}


##########
## main ##
##########

#test_uniform_law()
#test_normal_law()
#test_exponential_law()
#test_gamma_law()
#test_mixture_law()

## Parameters
loc = 0
scale = 1
shape = 1
size = 10000

## My implementation
rvX = RPOOstats::rv_GPD$new( loc = loc , scale = scale , shape = shape )

## Test rvs
if( FALSE )
{
	X = rvX$rvs(size)

	## Extremes packages
	Xe = extRemes::revd( n = size , loc = loc , scale = scale , shape = shape , threshold = loc , type = "GP" )
	
	p1 = ggplot2::ggplot( data.frame( X = X ) , ggplot2::aes( x = X ) ) + ggplot2::geom_histogram( ggplot2::aes( y = ..density.. ) , binwidth = 0.05 , color = "red" , fill = "red" , alpha = 0.5 ) + ggplot2::geom_density( alpha = 0.2 , color = "red" , fill = "red" ) + ggplot2::ggtitle( "Histogram" ) + ggplot2::xlim( base::c(0,10) )
	p2 = ggplot2::ggplot( data.frame( X = Xe ) , ggplot2::aes( x = X ) ) + ggplot2::geom_histogram( ggplot2::aes( y = ..density.. ) , binwidth = 0.05 , color = "red" , fill = "red" , alpha = 0.5 ) + ggplot2::geom_density( alpha = 0.2 , color = "red" , fill = "red" ) + ggplot2::ggtitle( "Histogram" )+ ggplot2::xlim( base::c(0,10) )
	gridExtra::grid.arrange( p1 , p2 , nrow = 1 )
}

## Test icdf
if( FALSE )
{
	q = seq( 0.01 , 0.99 , length = 100 )
	
	icdfX = rvX$icdf(q)
	icdfXe = qevd( q , loc = loc , scale = scale , shape = shape , threshold = loc , type = "GP" )
}

## Test cdf
if( FALSE )
{
	x = seq( 0 , 1 , length = 100 )
	cdfX = rvX$cdf(x)
	cdfXe = pevd( x , loc = loc , scale = scale , shape = shape , threshold = loc , type = "GP" )
}

## Test fit
if( TRUE )
{
	Xe = extRemes::revd( n = size , loc = loc , scale = scale , shape = shape , threshold = loc , type = "GP" )
	rvX = RPOOstats::rv_GPD$new()
	print( rvX$fit(Xe) )
}

