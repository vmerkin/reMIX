from numpy import arccos,tan,sin,cos
import point
import util

# theta prime angle
def theta_p(ps,pd):
        # this asserts that our calculation of theta_p agrees with
        # Vanhamaki, eq A3 assertion on floats is a bad idea anyway,
        # but I checked that the two calculations agree to the last
        # significant digit
        #
        # assert(cos(pd.t)*cos(ps.t)+sin(pd.t)*sin(ps.t)*cos(ps.p-pd.p)),
        # util.dot(ps.r_unit(),pd.r_unit())

	return arccos(util.dot(ps.r_unit(),pd.r_unit()))

# theta prime unit vector
def theta_p_unit(ps,pd):
	return pd.r_unit()/tan(theta_p(ps,pd))-ps.r_unit()/sin(theta_p(ps,pd))

# phi prime unit vector
def phi_p_unit(ps,pd):
	return util.cross(ps.r_unit(),pd.r_unit())/sin(theta_p(ps,pd))
