import ampere_pole,grid,matrices
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

    amp = ampere_pole.ampere_pole('data/AMPERE_north_2015-08-15T20-14-00Z.txt')

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

            A[K(i,jA),K(i,jA)] = -ft*(1/dt[i,j]+1/dt[i-1,j]) - fp*(1/dp[i,j]+1/dp[i,jm1])  - 0.5*ft*(dt[i,j]/dt[i-1,j]-dt[i-1,j]/dt[i,j])/tan(amp_grid.t[i,j])
            A[K(i,jA),K(i+1,jA)] = ft*(1./dt[i,j]+0.5*dt[i-1,j]/dt[i,j]/tan(amp_grid.t[i,j]))
            A[K(i,jA),K(i-1,jA)] = ft*(1./dt[i-1,j]-0.5*dt[i,j]/dt[i-1,j]/tan(amp_grid.t[i,j]))
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
        for jj in arange(0,Np): 
            A[K(0,j),K(1,jj)] = -dp[1,jj%(Np-1)]/(2.*pi)
        # note, not setting J because it's initializaed to zero anyway

    # Now, do the solve
    import scipy
    from scipy.sparse import linalg
    pot,status=scipy.sparse.linalg.lgmres(A,J)
    pot*=6.5**2*1.e3/Sigma

    print "Convergence status: ", status
    print "Potential min/max (kV): ",pot.min(),pot.max()

    # Plotting

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,44.*pi/180.)
    contourf(amp_grid.p+pi+pi/2,amp_grid.t,pot.reshape(24,51).T,51)
    colorbar().set_label('Potential, kV')

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,44.*pi/180.)
    contourf(amp_grid.p+pi+pi/2,amp_grid.t,amp_jr,51)
    colorbar().set_label('Current, microA/m^2')
    contour(amp_grid.p+pi+pi/2,amp_grid.t,pot.reshape(24,51).T,13,colors='black')
    

