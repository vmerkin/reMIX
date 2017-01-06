from numpy import sin,cos,arccos,mat,dot,repeat,multiply,tan,pi,sqrt,append,linspace,meshgrid,zeros,arange,ones_like,array,roll,hstack,newaxis,isscalar,ones,vstack
import scipy
from scipy import interpolate
from scipy.sparse import linalg
from scipy.sparse import coo_matrix

RIkm = 6.5e3  # Ionosphere radius in km

class solver():
    def __init__(self,grd):
        """ Initialize solver class.
        Arguments:
        grd -- grid object, see grid.py class
        """
        (self.Nt,self.Np) = grd.t.shape
        self.dt = grd.t[1:,:]-grd.t[:-1,:]
        self.dp = (roll(grd.p,-1,axis=1)-grd.p)%(2.*pi)    
        self.grd = grd
        # note, dp has the same size along phi axis is p: 
        # this eases the implementation of the periodic boundary below.

        self.cosd = self._cosDipAngle(grd.t)

    def _K(self,I,J):
        # Clean, elegant and pythonic
        return I[:,newaxis]*self.Np+J if not isscalar(I) else I*self.Np+J
        
    def _cosDipAngle(self,t):   # FIXME: modify to allow both hemispheres
        return(-2.*cos(t)/sqrt(1.+3.*cos(t)**2))

    def _defineInnerBlock(self,Nt,Np):
        return arange(1,Nt-1),arange(0,Np)

    def setLowLatBC(self,lowLatBC):
        self.lowLatBC = lowLatBC
    
    def setTerms(self,SigmaP):
        self.F11 = sin(self.grd.t)*SigmaP/self.cosd**2
        self.F22 = SigmaP
        
    def setStencilMatrixNp(self):
        # some aliases for brevity
        K = self._K
        Np = self.Np
        Nt = self.Nt
        dt = self.dt
        dp = self.dp
        F11 = self.F11
        F22 = self.F22
        grd = self.grd

        ############ inner block ############
        (I,J) = self._defineInnerBlock(Nt,Np)

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
        
        self.data = data
        self.I = II
        self.J = JJ

    def setRHSNp(self,S):
        """ Set the RHS for the matrix equation.
        
        Parameters:
        S -- the source term (i.e., radial current)
        """
        Np = self.Np
        Nt = self.Nt
        (I,J) = self._defineInnerBlock(Nt,Np)

        # Right hand side
        RHS = zeros(Nt*Np)
        RHS[self._K(I,J)] = S[I,:] # inner block
        RHS[self._K(Nt-1,J)] = self.lowLatBC # low lat BC
        # the rest are zeros

        self.RHS = RHS

    def setStencilMatrixAndRHS(self,S):
        # some aliases for brevity
        K = self._K
        Np = self.Np
        Nt = self.Nt
        dt = self.dt
        dp = self.dp
        F11 = self.F11
        F22 = self.F22
        grd = self.grd

        RHS = zeros(Nt*Np)
        nnz = Np*(Nt-2)*5 + Np + Np*(Np+1) # this is the number of non-zeros in the matrix. Note, this depends on the stencil and boundary conditions
        data=zeros(nnz)
        II = zeros(nnz)
        JJ = zeros(nnz)
        count = 0
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

                data[count] = -ft*( (F11[i,j]+F11[i+1,j])/dt[i,j]+(F11[i,j]+F11[i-1,j])/dt[i-1,j] ) - \
                              fp*( (F22[i,j]+F22[i,jp1])/dp[i,j%(Np-1)]+(F22[i,j]+F22[i,jm1])/dp[i, (j-1)%(Np-1)] )
                II[count] = K(i,j)
                JJ[count] = K(i,j)
                count+=1
                
                data[count] = ft*(F11[i,j]+F11[i+1,j])/dt[i,j]
                II[count]  = K(i,j)
                JJ[count] = K(i+1,j)
                count+=1

                data[count] = ft*(F11[i,j]+F11[i-1,j])/dt[i-1,j]
                II [count] = K(i,j)
                JJ[count]  = K(i-1,j)
                count+=1

                data[count] = fp*(F22[i,j]+F22[i,jp1])/dp[i, j%(Np-1)]
                II[count] = K(i,j)
                JJ[count] = K(i,jp1)
                count+=1

                data[count] = fp*(F22[i,j]+F22[i,jm1])/dp[i, (j-1)%(Np-1)]
                II[count] = K(i,j)
                JJ[count] = K(i,jm1)
                count+=1

                RHS[K(i,j)] = S[i,j]

        # low lat boundary
        for j in arange(0,Np):
            data[count]  = 1.
            II[count] = K(Nt-1,j)
            JJ[count] = K(Nt-1,j)
            count+=1

            RHS[K(Nt-1,j)] = 0. # FIXME: make arbitrary BC

        # pole boundary
        for j in arange(0,Np):
            data[count] = 1.
            II[count] = K(0,j)
            JJ[count] =  K(0,j)
            count+=1

            for jj in arange(0,Np): 
                data[count] = -dp[1,jj%(Np-1)]/(2.*pi)
                II[count] = K(0,j)
                JJ[count] = K(1,jj)
                count+=1
                # note, not setting J because it's initializaed to zero anyway

        self.data = data
        self.I = II
        self.J = JJ
        self.RHS = RHS

    def solve(self):
        Np = self.Np
        Nt = self.Nt
        
        M = coo_matrix( (self.data,(self.I,self.J)),shape=(Np*Nt,Np*Nt) )
        #    pot,status=scipy.sparse.linalg.lgmres(M,Jr)    # this becomes very slow for large matrices, need to precondition?
        #    print "Convergence status: ", status

        #    solve=scipy.sparse.linalg.factorized(M)   # this is similar to above
        #    pot=solve(Jr)

        pot=scipy.sparse.linalg.spsolve(M.tocsc(),self.RHS)
        return (pot*RIkm**2*1.e-3)  # assume current was in microA/m^2


    
