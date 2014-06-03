########################################################################
# R.A.Borrelli
# @TheDoctorRAB
# rev.31.May.2014
########################################################################
# FEM solution to the diffusion reaction equation.
# These are the required functions.
# These are indexed by letters and referenced in the main program.
# Because it is just easy to search for them if needed.
# The names may seem overly much.
# I used them in my dissertation code which was in FORTRAN.
########################################################################
#
########################################################################
# e,i,j,k,n are used for loop indices
# e=element based loops
# i,j,k=nodal based loops
# n=time step based loops
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
####### (a): evaluate linear function
# y=m*x+b
# m=slope
# x=position in real space (x) or isoparametric space (xi)
# b=intercept
# y=result
# x is used here in the function for convenience
###
#
###
def evaluation_linear(m,x,b):
###
    y=m*x+b
    return y
#######
#
#
#
####### (b): Gaussian quadrature for numerical integration
# Numerical integration is performed by applying gaussian quadrature
# This is for linear functions of form m*x+b in the isoparametric space
#
# So the linear function to be integrated must be transposed to the isoparametric space first
#
#             *******ARGUMENT ORDER IS CRITICAL*******
#
# 01: jacobian=jacobian to map from real space to iso-space
# 02: slope_m1=slope (m) for shape function 1
# 03: intercept_b1=intercept (b) for shape function 1
# 04: slope_m2=slope (m) for shape function 2
# 05: intercept_b2=intercept (b) for shape function 2
# 06: function_slope_m=slope (m) for function f(x)
# 07: function_intercept_b=intercept (b) for function f(x)
# 08: isomap_slope_m=slope (m) for the isoparametric function x(xi)
# 09: isomap_intercept_b=intercept (b) for the isoparametric function x(xi)
#
# if there is a deriviative of the shape function needed then put the derivative in the b slot (i.e., +/- 1/2) and m=0
# if there is no function f(x), then put m=0,b=1 in 08,09
###
#
###
def gaussian_quadrature_linear(jacobian,slope_m1,intercept_b1,slope_m2,intercept_b2,function_slope_m,function_intercept_b,isomap_slope_m,isomap_intercept_b):
###
# read in the gaussian quadrature points
# the input file should be structured [gauss_points,2]
    quadrature=numpy.loadtxt('quadrature.inp')
###
# determine the number of gauss points    
# initialize counter to count the number of students
    i=0
###
# loop through to count the students
    for line in quadrature:
        i=i+1
        gauss_points=i
# end line
###
# initialize the addition storage variable
    integral_response=float(0)
###
# initialize local quadrature
    local_quadrature=numpy.zeros((gauss_points))
###
# perform numerical integration
    for j in range(0,gauss_points):
        local_quadrature[j]=quadrature[j,0]*(2*jacobian*evaluation_linear(slope_m1,intercept_b1,quadrature[j,1])*evaluation_linear(slope_m2,intercept_b2,quadrature[j,1])*evaluation_linear(function_slope_m,function_intercept_b,quadrature[j,1])*evaluation_linear(isomap_slope_m,isomap_intercept_b,quadrature[j,1]))
# end j
###
# add the integration
    for k in range(0,gauss_points):
        integral_response=integral_response+local_quadrature[k]
# end k
###
    return integral_response
#######
#
#
#
####### (c): spatial mesh dicretization
# This is a uniform mesh over the entire time domain
# The spatial positions of the nodes and element lengths are computed
# TRIBEX has the code for temporal dependent mesh with different resolutions
###
#
###
def spatial_mesh_discretization():
###
# read in mesh parameters
    mesh_parameters=numpy.loadtxt('fem_diffusion.equation_mesh.parameters.inp')
###
# assign spatial parameters
    x0=mesh_parameters[0]
    xSink=mesh_parameters[1]
###
# assign mesh parameters
    node=int(mesh_parameters[2])
    element=node-1
###
# initialize the spatial mesh and element lengths
    spatial_mesh=numpy.zeros((node))
    h=numpy.zeros((element))
###
# build the spatial mesh
# node 1
    spatial_mesh[0]=x0
###
# and the rest
    for i in range(1,node):
        spatial_mesh[i]=spatial_mesh[i-1]+(xSink-x0)/(element)
# end i
###
# calculate the element lengths
    for e in range(0,element):
        h[e]=spatial_mesh[e+1]-spatial_mesh[e]
# end e  
###
    return (spatial_mesh,h,node,element)
#######
#
#
#
########################################################################
#      EOF
########################################################################
