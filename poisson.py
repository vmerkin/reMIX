import ampere_pole,grid,matrices
from numpy import linalg as la
import numpy as np
import sys
import time
import point,sphtrig
from numpy import sin,cos,arccos,mat,dot,repeat,multiply,tan,pi,sqrt,append,linspace,meshgrid,zeros,arange,ones_like,array
from matplotlib.pyplot import figure,subplot,pcolormesh,colorbar,contour,show,contourf,rgrids,ylim
from scipy import interpolate
import time

if __name__ == "__main__":
    Sigma_const = 8. # scalar for the time being
    interp = False

    # AMPERE stuff
    amp = ampere_pole.ampere_pole('data/AMPERE_north_2015-08-15T20-50-00Z.txt')

    grd = grid.generate_from_data(amp.get_grid())
    jr = amp.get_data()
    if interp:
        ntheta = 100
        nphi = 720
        f = interpolate.RectBivariateSpline(grd.t[:,0],grd.p[0,:],jr,kx=1,ky=1)
        new_t = linspace(grd.t[:,0].min(),grd.t[:,0].max(),ntheta)
        new_p = linspace(grd.p[0,:].min(),grd.p[0,:].max(),nphi)
        grd = grid.grid()
        grd.t,grd.p = meshgrid(new_t,new_p,indexing='ij')
        jr = f(grd.t[:,0],grd.p[0,:])

    jr[abs(jr)<=0.2]=0.

    dt = grd.t[1:,:]-grd.t[:-1,:]
    dp = grd.p[:,1:]-grd.p[:,:-1]

    (Nt,Np) = grd.t.shape
    ############

    # define conductance as uniform Pedersen for now
    SigmaP = ones_like(jr)*Sigma_const
    cosd     = -2.*cos(grd.t)/sqrt(1.+3.*cos(grd.t)**2)
    F11          = sin(grd.t)*SigmaP/cosd**2
    F22          = SigmaP


    K = lambda i,j: j*Nt+i
    A=zeros((Nt*Np,Nt*Np))
    J = zeros(Nt*Np)

    # inner block
    for i in arange(1,Nt-1):
        for j in arange(0,Np):
            jm1 = (j-1)%Np
            jp1  = (j+1)%Np
            j     = j%Np  # not needed (J<=Np-1) but do it for symmetry

            # note, we use the above definitions of j's except in reference to dp array, which has Np-1 points
            # in those cases, we take the modulo in place

            ft = 1./(dt[i,j]+dt[i-1,j])/sin(grd.t[i,j])
            fp = 1./(dp[i,j%(Np-1)]+dp[i, (j-1)%(Np-1)])/sin(grd.t[i,j])**2

            A[K(i,j),K(i,j)] = -ft*( (F11[i,j]+F11[i+1,j])/dt[i,j]+(F11[i,j]+F11[i-1,j])/dt[i-1,j] ) - \
                                 fp*( (F22[i,j]+F22[i,jp1])/dp[i,j%(Np-1)]+(F22[i,j]+F22[i,jm1])/dp[i, (j-1)%(Np-1)] )
            A[K(i,j),K(i+1,j)] = ft*(F11[i,j]+F11[i+1,j])/dt[i,j]
            A[K(i,j),K(i-1,j)] = ft*(F11[i,j]+F11[i-1,j])/dt[i-1,j]
            A[K(i,j),K(i,jp1)] = fp*(F22[i,j]+F22[i,jp1])/dp[i, j%(Np-1)]
            A[K(i,j),K(i,jm1)] = fp*(F22[i,j]+F22[i,jm1])/dp[i, (j-1)%(Np-1)]

            J[K(i,j)] = jr[i,j]


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
    pot*=6.5**2*1.e3

    print "Convergence status: ", status
    print "Potential min/max (kV): ",pot.min(),pot.max()

    # Plotting

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,44.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,pot.reshape(Np,Nt).T,51)
    colorbar().set_label('Potential, kV')

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,44.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,jr,51)
    colorbar().set_label('Current, microA/m^2')
    contour(grd.p+pi+pi/2,grd.t,pot.reshape(Np,Nt).T,13,colors='black')
    

