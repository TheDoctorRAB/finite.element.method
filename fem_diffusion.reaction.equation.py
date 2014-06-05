########################################################################
# R.A.Borrelli
# @TheDoctorRAB
# rev.05.June.2014
# v1.0.0
########################################################################
# FEM solution to the diffusion reaction equation (time-dependent).
# Transport in a porous medium.
# Constant porosity, diffusion coefficient, sorption distribution.
# Cartesian coordinate system.
# Semi-infinite domain.
# Dirichlet boundary conditions.
# Initial condition is zero.
########################################################################
#
########################################################################
# The mesh is currently uniform and the same at all time steps.
# This can be improved on with the TRIBEX moving mesh later.
########################################################################
#
########################################################################
# The processor is a preconditioned conjugate gradient (PCG).
# Made from scratch based on the TRIBEX model.
########################################################################
#
########################################################################
# The first two time steps used a backward Euler method.
# Then an implicit trapezoid methos is used.
########################################################################
#
########################################################################
# Loop convention
###
# e,i,j,k,n are used for loop indices
# e=element based loops
# i,j,k=nodal based loops
# n=time step based loops
########################################################################
#
########################################################################
# Domain
###
# spatial_mesh=spatial domain
# h=length of each element in the spatial domain
# node=total number of nodes
# element=total number of elements
########################################################################
#
########################################################################
# Time
###
# tdomain=total number of time steps ***includes t=0***
# full_step=time domain less t=0, so full_step=tdomain-1
########################################################################
#
########################################################################
# Physical domain
###
# x0=origin; left boundary
# xSink=sink; right boundary
# xSink has to be 'large' to describe a semi-infinite domain
########################################################################
#
########################################################################
# Environmental parameters
###
# half_life=radionuclide half-life
# diffusion_coefficient=molecular diffusion coefficient in free water
# bulk_density=bentonite bulk density
# sorption_distribution=sorption distribution coefficient
# porosity=fluid porosity
# decay constant=decay constant; reaction constant
# sorption=dimensionaless constant used in retardaction factor
# solid_volume=solid volume fraction
# retardation=retardation factor
# bc1,bc2=boundary conditions at x=0,Sink
########################################################################
#
########################################################################
# FEM
###
# local_parametric_mapping=transfroms real space position to iso-space
# shape_function=linear functions in iso-space
# dshape_functions=derivatives of shape functions
###
# threshold,delta,chi,epsilon,numerator,denominator=PCG procedurals
########################################################################
#
########################################################################
# Gaussian integration
###
# gauss_points=number of points in the integration
# quaradture=integration points and functional values
########################################################################
#
#
#
####### imports
import numpy
import side_conditions
import fem_functions as fem_f
import physical_constants as phyc
######
#
#
#
########################################################################
#
#
#
###### preprocessing #######
#
###### read in the input data
# boundary conditions
bc1,bc2=side_conditions.side_conditions()
###
# physical constants
decay_constant,diffusion_coefficient,porosity,retardation=phyc.physical_constants()
###
# time steps
tdomain,full_step,time_domain=fem_f.time_steps()
###
# construct the spatial domain
spatial_mesh,h,node,element=fem_f.spatial_mesh_discretization()
###
# initialize the solution matrix
radionuclide_concentration=fem_f.initialize_solution_matrix(full_step,node,spatial_mesh,time_domain)
numpy.savetxt('sol.out',radionuclide_concentration)

























