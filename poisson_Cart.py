import ampere,grid,matrices
from numpy import linalg as la
import numpy as np
import sys
import time
import point,sphtrig
from numpy import sin,cos,arccos,mat,dot,repeat,multiply,tan,pi,sqrt,append,linspace,meshgrid,zeros,arange
from matplotlib.pyplot import figure,subplot,pcolormesh,colorbar,contour,show,contourf,rgrids,ylim
from scipy import interpolate

if __name__ == "__main__":
    Sigma = 8. # scalar for the time being

    amp = ampere.ampere('data/AMPERE_north_2015-08-15T12-00-00Z.txt')

    amp_grid = grid.generate_from_data(amp.get_grid())
    amp_jr = amp.get_data()

    dt = amp_grid.t[1:,:]-amp_grid.t[:-1,:]
    dp = amp_grid.p[:,1:]-amp_grid.p[:,:-1]

    (Nt,Np) = amp_grid.t.shape

    K = lambda i,j: j*Nt+i
    A=zeros((Nt*Np,Nt*Np))
    J = zeros(Nt*Np)

    # inner block
    for i in arange(1,Nt-1):
        for j in arange(0,Np):
            jAm1 = (j-1)%Np
            jAp1 = (j+1)%Np
            jA = j

            jm1 = (j-1)%(Np-1)
            jp1  = (j+1)%(Np-1)
            j      = j%(Np-1)

            ft = 2./(dt[i,j]+dt[i-1,j])
            fp = 2./(dp[i,j]+dp[i,jm1])/sin(amp_grid.t[i,j])**2

            A[K(i,jA),K(i,jA)] = -ft*(1/dt[i,j]+1/dt[i-1,j]) - fp*(1/dp[i,j]+1/dp[i,jm1])
            A[K(i,jA),K(i+1,jA)] = ft/dt[i,j]
            A[K(i,jA),K(i-1,jA)] = ft/dt[i-1,j]
            A[K(i,jA),K(i,jAp1)] = fp/dp[i,j]
            A[K(i,jA),K(i,jAm1)] = fp/dp[i,jm1]

            J[K(i,jA)] = amp_jr[i,j]

    # low lat boundary
    for j in arange(0,Np):
        A[K(Nt-1,j),K(Nt-1,j)] = 1.
        J[K(Nt-1,j)] = 0. # FIXME: make arbitrary BC

    # pole boundary
    for j in arange(0,Np):
        A[K(0,j),K(0,j)] = 1.
        J[K(0,j)] = 0. # FIXME: enforce proper BC




    # Plotting

    # fixup periodic boundary before plotting
    



    # e_plot = append(e,e[:,[0]],axis=1)
    # phi_plot = append(phi,phi[:,[0]],axis=1)
    # p_plot = append(g1.p,g1.p[:,[0]]+2*pi,axis=1)
    # t_plot = append(g1.t,g1.t[:,[0]],axis=1)
    # t1_plot = append(amp_grid.t,amp_grid.t[:,[0]],axis=1)


    # figure()
    # subplot(111,polar=True)
    # circles = linspace(10,40,4)*pi/180.
    # lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    # lines, labels = rgrids(circles,lbls)
    # ylim(0,44.*pi/180.)
    # contourf(amp_grid.p+pi+pi/2,amp_grid.t,amp_jr,51)
#    colorbar().set_label('|E|, mV/m')
