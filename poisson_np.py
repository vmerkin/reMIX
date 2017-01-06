import ampere_pole,grid,matrices
from numpy import linalg as la
import numpy as np
import sys
import time
import point,sphtrig
from numpy import sin,cos,arccos,mat,dot,repeat,multiply,tan,pi,sqrt,append,linspace,meshgrid,zeros,arange,ones_like,array,roll,hstack,newaxis,isscalar,ones,vstack
from matplotlib.pyplot import figure,subplot,pcolormesh,colorbar,contour,show,contourf,rgrids,ylim
import scipy
from scipy import interpolate
from scipy.sparse import linalg
from scipy.sparse import coo_matrix
import time

if __name__ == "__main__":
    Sigma_const = 8. # scalar for the time being
    interp = True

    # AMPERE stuff
    amp = ampere_pole.ampere_pole('data/AMPERE_north_2015-08-15T20-50-00Z.txt')

    grd = grid.generate_from_data(amp.get_grid())
    jr = amp.get_data()
    if interp:
        ntheta = 50
        nphi = 360

        # Complete full circle with AMPERE data so that interpolated grid doesn't have a huge hole
        # note, however, that we still perform the solve on a cell-centered grid, i.e., in the case of interpolation
        # here we make sure that the resultant grid doesn't wrap around
        p_tmp = hstack((grd.p,2.*pi+grd.p[:,[0]]))
        t_tmp = hstack((grd.t,grd.t[:,[0]]))
        jr_tmp = hstack((jr,jr[:,[0]]))

        f = interpolate.RectBivariateSpline(grd.t[:,0],grd.p[0,:],jr,kx=1,ky=1)  # note, assume uniform ampere grid here
        # this here defines cell vortices (and wraps around)
        new_t = linspace(t_tmp.min(),t_tmp.max(),ntheta)
        new_p = linspace(p_tmp.min(),p_tmp[0,:].max(),nphi+1)
        # now redefine as cell centers
        new_p = 0.5*(new_p[:-1]+new_p[1:])
        grd = grid.grid()
        grd.t,grd.p = meshgrid(new_t,new_p,indexing='ij')
        jr = f(grd.t[:,0],grd.p[0,:])

    jr[abs(jr)<=0.2]=0.

    dt = grd.t[1:,:]-grd.t[:-1,:]
    dp = (roll(grd.p,-1,axis=1)-grd.p)%(2.*pi)

    (Nt,Np) = grd.t.shape
    ############

    # define conductance as uniform Pedersen for now
    SigmaP = ones_like(jr)*Sigma_const

    # low lat BC values, set to zero for now
    lowLatBC = zeros(Np)

    cosd     = -2.*cos(grd.t)/sqrt(1.+3.*cos(grd.t)**2)
    F11          = sin(grd.t)*SigmaP/cosd**2
    F22          = SigmaP

    t0=time.clock()

    # Clean, elegant and pythonic
    K = lambda I,J: I[:,newaxis]*Np+J if not isscalar(I) else I*Np+J

    ############ inner block ############
    I = arange(1,Nt-1)
    J = arange(0,Np)

    # This works but not elegant. Use the above definition of K lambda, but saving it here for the record
    # IJ = meshgrid(I,J,indexing='ij')
    # Kij   = IJ[0]*Np+IJ[1]%Np

    ft = 1./(dt[I,:]+dt[I-1,:])/sin(grd.t[I,:])
    fp = 1./(dp[I,:]+roll(dp[I,:],1,axis=1))/sin(grd.t[I,:])**2

    ijij   = -ft*( (F11[I,:]+F11[I+1,:])/dt[I,:]+(F11[I,:]+F11[I-1,:])/dt[I-1,:] ) - \
             fp*( (F22[I,:]+roll(F22[I,:],-1,axis=1))/dp[I,:]+(F22[I,:]+roll(F22[I,:],1,axis=1))/roll(dp[I,:],1,axis=1) )
    ijip1j = ft*(F11[I,:]+F11[I+1,:])/dt[I,:] 
    ijim1j = ft*(F11[I,:]+F11[I-1,:])/dt[I-1,:]
    ijijp1 = fp*(F22[I,:]+roll(F22[I,:],-1,axis=1))/dp[I,:]
    ijijm1 = fp*(F22[I,:]+roll(F22[I,:],1,axis=1))/roll(dp[I,:],1,axis=1) 

    Kij   = K(I,J)
    Kip1j = K(I+1,J)
    Kim1j = K(I-1,J)
    Kijp1 = K(I,(J+1)%Np)
    Kijm1 = K(I,(J-1)%Np)
    ############ inner block ############

    data = hstack(
        (ijij.ravel(),
         ijip1j.ravel(),
         ijim1j.ravel(),
         ijijp1.ravel(),
         ijijm1.ravel(),
         ones(Np),    # low lat boundary
         ones(Np),    # pole boundary
         hstack(Np*[-dp[1,:]/(2.*pi)]))
    )
    
    II = hstack(
        (Kij.ravel(),
         Kij.ravel(),
         Kij.ravel(),
         Kij.ravel(),
         Kij.ravel(),
         K(Nt-1,J),  # low lat boundary
         K(0,J),     # pole boundary
         K(zeros(Np),J).T.ravel())
#         hstack([ones(Np)*K(0,j) for j in J]))  # this is slightly slower than the above, but the results are identical
    )

    JJ = hstack(
        (Kij.ravel(),
         Kip1j.ravel(),
         Kim1j.ravel(),
         Kijp1.ravel(),
         Kijm1.ravel(),
         K(Nt-1,J), # low lat boundary
         K(0,J),    # pole boundary
         hstack(Np*[K(1,J)]))
        )

    # Right hand side
    Jr = zeros(Nt*Np)
    Jr[K(I,J)] = jr[I,:] # inner block
    Jr[K(Nt-1,J)] = lowLatBC # low lat BC
    # the rest are zeros

    print('Time spend constructing matrix (s):',time.clock()-t0)
    t0 = time.clock()
    # Now, do the solve
    M = coo_matrix( (data,(II,JJ)),shape=(Np*Nt,Np*Nt) )
#    pot,status=scipy.sparse.linalg.lgmres(M,Jr)    # this becomes very slow for large matrices, need to precondition?
#    solve=scipy.sparse.linalg.factorized(M)   # this is similar to above
#    pot=solve(Jr)

    pot=scipy.sparse.linalg.spsolve(M.tocsc(),Jr)

    print('Time spend on solve (s):',time.clock()-t0)

    pot*=6.5**2*1.e3

#    print "Convergence status: ", status
    print "Potential min/max (kV): ",pot.min(),pot.max()

    # Plotting

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,50.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,pot.reshape(Nt,Np),51)
    colorbar().set_label('Potential, kV')

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,50.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,jr,45)
    colorbar().set_label('Current, microA/m^2')
    contour(grd.p+pi+pi/2,grd.t,pot.reshape(Nt,Np),13,colors='black')
    

