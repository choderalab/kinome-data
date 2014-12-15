#!/usr/bin/env python

"""
A test of pymc for ITC.

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy
import copy

import scipy.optimize
import scipy.integrate

from math import sqrt, exp, log

#=============================================================================================
# Physical constants
#=============================================================================================

Na = 6.02214179e23 # Avogadro's number (number/mol)
kB = Na * 1.3806504e-23 / 4184.0 # Boltzmann constant (kcal/mol/K)
C0 = 1.0 # standard concentration (M)

#=============================================================================================
# Competitive binding model
#=============================================================================================

class CompetitiveBindingModel(object):
   """
   Competitive binding model.

   """

   @classmethod
   def equilibrium_concentrations(cls, Ka_n, C0_R, C0_Ln, V, c0=None):
      """
      Compute the equilibrium concentrations of each complex species for N ligands competitively binding to a receptor.

      ARGUMENTS

      Ka_n (numpy N-array of float) - Ka_n[n] is the association constant for receptor and ligand species n (1/M)
      x_R (float) - the total number of moles of receptor in the sample volume
      x_n (numpy N-array of float) - x_n[n] is the total number of moles of ligand species n in the sample volume
      V (float) - the total sample volume (L)

      RETURNS

      C_n (numpy N-array of float) - C_n[n] is the concentration of complex of receptor with ligand species n

      EXAMPLES

      >>> V = 1.4303e-3 # volume (L)
      >>> x_R = V * 510.e-3 # receptor
      >>> x_Ln = numpy.array([V * 8.6e-6, 200.e-6 * 55.e-6]) # ligands
      >>> Ka_n = numpy.array([1./(400.e-9), 1./(2.e-11)]) # association constants
      >>> C_PLn = equilibrium_concentrations(Ka_n, x_R, x_Ln, V)

      NOTES
      
      Each complex concentration C_n must obey the relation
      
      Ka_n[n] = C_RLn[n] / (C_R * C_Ln[n])           for n = 1..N
      
      with conservation of mass constraints
      
      V * (C_Ln[n] + C_RLn[n]) = x_Ln[n]             for n = 1..N


      
      and
      
      V * (C_R + C_RLn[:].sum()) = x_R
      
      along with the constraints
      
      0 <= V * C_RLn[n] <= min(x_Ln[n], x_R)         for n = 1..N
      V * C_RLn[:].sum() <= x_R
      
      We can rearrange these expressions to give
      
      V * C_R * C_Ln[n] * Ka_n[n] - V * C_RLn[n] = 0
      
      and eliminate C_Ln[n] and C_R to give
      
      V * (x_R/V - C_RLn[:].sum()) * (x_Ln[n]/V - C_RLn[n]) * Ka_n[n] - V * C_RLn[n] = 0    for n = 1..N
      
      """

      x_R = C0_R * V
      x_Ln = C0_Ln * V
      
      nspecies = Ka_n.size
      #print "x_R = ", x_R
      #print "x_Ln = ", x_Ln
      #print "x_Ln / V = ", x_Ln / V
      #print "Ka_n = ", Ka_n

      # Define optimization functions
      def func(C_RLn):
         f_n = V * (x_R/V - C_RLn[:].sum()) * (x_Ln[:]/V - C_RLn[:]) * Ka_n[:] - V * C_RLn[:]
         #print "f_n = ", f_n
         return f_n
      
      def fprime(C_RLn):
         nspecies = C_RLn.size
         G_nm = numpy.zeros([nspecies,nspecies], numpy.float64) # G_nm[n,m] is the derivative of func[n] with respect to C_RLn[m]
         for n in range(nspecies):
            G_nm[n,:] = - V * (x_Ln[:]/V - C_RLn[:]) * Ka_n[:]
            G_nm[n,n] -= V * (Ka_n[n] * (x_R/V - C_RLn[:].sum()) + 1.0)
         return G_nm

      def sfunc(s):
         #print "s = ", s
         f_n = V * (x_R/V - (s[:]**2).sum()) * (x_Ln[:]/V - s[:]**2) * Ka_n[:] - V * s[:]**2
         #print "f_n = ", f_n
         return f_n
      
      def sfprime(s):
         nspecies = s.size
         G_nm = numpy.zeros([nspecies,nspecies], numpy.float64) # G_nm[n,m] is the derivative of func[n] with respect to C_RLn[m]
         for n in range(nspecies):
            G_nm[n,:] = - V * (x_Ln[:]/V - s[:]**2) * Ka_n[:] 
            G_nm[n,n] -= V * (Ka_n[n] * (x_R/V - (s[:]**2).sum()) + 1.0)
            G_nm[n,:] *= 2. * s[n]
         return G_nm

      # Allocate storage for complexes
      # Compute equilibrium concentrations.
      #x0 = numpy.zeros([nspecies], numpy.float64)
      #x0 = (x_Ln / V).copy()
      #x = scipy.optimize.fsolve(func, x0, fprime=fprime)
      #C_RLn = x
      
      #x0 = numpy.sqrt(x_Ln / V).copy()
      #x = scipy.optimize.fsolve(sfunc, x0, fprime=sfprime)         
      #C_RLn = x**2

      def objective(x):
         f_n = func(x)
         G_nm = fprime(x)

         obj = (f_n**2).sum()
         grad = 0.0 * f_n
         for n in range(f_n.size):
            grad += 2 * f_n[n] * G_nm[n,:]

         return (obj, grad)
         
      #x0 = numpy.zeros([nspecies], numpy.float64)      
      #bounds = list()
      #for n in range(nspecies):
      #   m = min(C0_R, C0_Ln[n])
      #   bounds.append( (0., m) )
      #[x, a, b] = scipy.optimize.fmin_l_bfgs_b(objective, x0, bounds=bounds)
      #C_RLn = x

      def ode(c_n, t, Ka_n, x_Ln, x_R):
         dc_n = - c_n[:] + Ka_n[:] * (x_Ln[:]/V - c_n[:]) * (x_R/V - c_n[:].sum())
         return dc_n
      
      def odegrad(c_n, t, Ka_n, x_Ln, x_R):
         N = c_n.size
         d2c = numpy.zeros([N,N], numpy.float64)
         for n in range(N):
            d2c[n,:] = -Ka_n[n] * (x_Ln[n]/V - c_n[n])
            d2c[n,n] += -(Ka_n[n] * (x_R/V - c_n[:].sum()) + 1.0)
         return d2c

      #if c0 is None: c0 = numpy.zeros([nspecies], numpy.float64)
      #maxtime = 100.0 * (x_R/V) / Ka_n.max()
      #time = [0, maxtime / 2.0, maxtime]
      #c = scipy.integrate.odeint(ode, c0, time, Dfun=odegrad, args=(Ka_n, x_Ln, x_R))
      #C_RLn = c[-1,:]

      #c = numpy.zeros([nspecies], numpy.float64)
      #maxtime = 1.0 / Ka_n.min()
      #maxtime = 1.0 / ((x_R/V) * Ka_n.min())
      #maxtime = 1.0
      #time = [0, maxtime]
      #c = scipy.optimize.fsolve(ode, c, fprime=odegrad, args=(0.0, Ka_n, x_Ln, x_R), xtol=1.0e-6)
      #c = scipy.integrate.odeint(ode, c, time, Dfun=odegrad, args=(Ka_n, x_Ln, x_R), mxstep=50000)
      #c = c[-1,:]
      #C_RLn = c

      #print "C_RLn = ", C_RLn 
      #print ""

      c = numpy.zeros([nspecies], numpy.float64)
      sorted_indices = numpy.argsort(-x_Ln)
      for n in range(nspecies):
         indices = sorted_indices[0:n+1]
         #c[indices] = scipy.optimize.fsolve(ode, c[indices], fprime=odegrad, args=(0.0, Ka_n[indices], x_Ln[indices], x_R), xtol=1.0e-6, warning=False)
         c[indices] = scipy.optimize.fsolve(ode, c[indices], fprime=odegrad, args=(0.0, Ka_n[indices], x_Ln[indices], x_R), xtol=1.0e-6)
      C_RLn = c
      
      return C_RLn
      
