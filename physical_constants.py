########################################################################
# R.A.Borrelli
# @TheDoctorRAB
# rev.02.June.2014
########################################################################
# This reads in physicial constants needed for the diffusion equation.
# The file is two columns: name,value.
########################################################################
# e,i,j,k,n are used for loop indices
########################################################################
#
#
#
####### imports
import numpy
#######
#
#
#
####### read and compute physical constants
###
def physical_constants():
###
# read in the data
    physical_constants=numpy.loadtxt('physical_constants.inp',dtype=str)
    half_life=float(physical_constants[0,1])
    diffusion_coefficient=float(physical_constants[1,1])
    bulk_density=float(physical_constants[2,1])
    sorption_distribution=float(physical_constants[3,1])
    porosity=float(physical_constants[4,1])
###
# compute the decay constant in 1/s
    decay_constant=(numpy.log(2))/(half_life*365.25*24*3600)
###
# compute sorption (K in the TRIBEX model)
    sorption=bulk_density*sorption_distribution
###
# compute solid_volume_fraction
    solid_volume_fraction=1-porosity
###
# compute 'retardation factor'
# this is not the formal definition of retardation factor
# but actually fluid porosity*retardation factor
# there is not any name for it
    retardation=porosity+solid_volume_fraction*sorption
###
    return (decay_constant,diffusion_coefficient,porosity,retardation)
#######
#
#
#
########################################################################
#      EOF
########################################################################
