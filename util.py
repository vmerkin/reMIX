import numpy as np

# dot product of two 3D vectors
def dot(v1,v2):
	assert v1.size==3 and v2.size==3, 'Both vectors should have size=3'

	return(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

# cross product of two 3D vectors
def cross(v1,v2):
	assert v1.size==3 and v2.size==3, 'Both vectors should have size=3'

	return(np.array([v1[1]*v2[2]-v1[2]*v2[1],
                      v1[2]*v2[0]-v1[0]*v2[2],
                      v1[0]*v2[1]-v1[1]*v2[0]]))
        
# magnitude of a vector
def mag(v):
	assert v.size==3 , 'Vector should have size=3'

	return(np.sqrt(v[0]**2+v[1]**2+v[2]**2))
