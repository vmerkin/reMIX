from numpy import arccos,tan,sin
import point
import util

# theta prime angle
def t_p(ps,po):
	return arccos(util.dot(ps.r_u(),po.r_u()))

# theta prime unit vector
def t_p_u(ps,po):
	return po.r_u()/tan(t_p(ps,po))-ps.r_u()/sin(t_p(ps,po))

# phi prime unit vector
def p_p_u(ps,po):
	return util.cross(ps.r_u(),po.r_u())/sin(t_p(ps,po))
