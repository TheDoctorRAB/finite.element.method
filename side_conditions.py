########################################################################
# R.A.Borrelli
# @TheDoctorRAB
# rev.03.June.2014
########################################################################
# Reads in boundary conditions and time for the diffusion reaction equation.
# side_conditions: condition at x=0, condition at x=sink, time (s)
# Enter one per line.
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
####### read in the conditions
def side_conditions():
    conditions=numpy.loadtxt('side_conditions.inp')
    bc1=conditions[0]
    bc2=conditions[1]
    time=conditions[2]
    return (bc1,bc2,time)
#######
#
#
#
########################################################################
#      EOF
########################################################################
