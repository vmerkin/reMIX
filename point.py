from numpy import cos,sin,array,pi

class point():

    def __init__(self,t,p):
    	# assume both theta and phi are in radians
    	assert t<=pi/2 and t>=0. and p<=pi and p>=0.
        self.t  = t
        self.p = p

  	# unit vector
    def r_u(self):
    	return(array([sin(self.t)*cos(self.p),sin(self.t)*sin(self.p),cos(self.t)]))