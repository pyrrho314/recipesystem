import numpy as np
import math

def rms(points1, points2, axis=None):
    if axis is not None:
        points1 = points1[:,axis]
        points2 = points2[:,axis]

    diff = points1 - points2
    npoints = len(points1) - 1
    return np.sqrt(np.sum(diff**2) / npoints)

def piecewise_rms(left, rght):
    x_rms = rms(left, rght, 0)
    y_rms = rms(left, rght, 1)
    return np.array([x_rms, y_rms])

def rmsd(rms):
    return np.sqrt(np.sum(rms**2))

class Transformation(object):
    def __init__(self, euler_angle=None, rotmat=None, translation=None):
        self.matrix = np.identity(3)

        if euler_angle is not None:
            self.push_front_euler_rotation(euler_angle)
        if rotmat is not None:
            self.push_front_rotation(rotmat)
        if translation is not None:
            self.push_front_translation(translation)

    def push_front(self, transform):
        if isinstance(transform, np.ndarray):
            assert(transform.shape == (3, 3))
            self.matrix = transform.dot(self.matrix)
        else:
            self.push_front(transform.matrix)

    def push_front_translation(self, trans):
        assert(trans.shape == (2,))
        matrix = np.identity(3)
        matrix[:2,-1] = trans
        self.push_front(matrix)

    def push_front_rotation(self, rmat):
        assert(rmat.shape == (2, 2))
        matrix = np.identity(3)
        matrix[:2,:2] = rmat
        self.push_front(matrix)

    def push_front_euler_rotation(self, deg):
        """ Counter-clockwise rotation angle in degrees """
        rad = math.radians(deg)
        rmat = np.array([[math.cos(rad), -math.sin(rad)],
                         [math.sin(rad),  math.cos(rad)]])
        self.push_front_rotation(rmat)

    def transform(self, points):
        assert(self.matrix.shape == (3, 3))
        assert(points.shape[1] == 2)

        points3d = np.ones((points.shape[0], 3))
        points3d[:,:2] = points
        transformed = np.dot(points3d, self.matrix.transpose())

        return transformed[:,:2]

    def get_euler_angle(self):
        assert(self.matrix.shape == (3, 3))
        return math.degrees(math.atan2(-self.matrix[0][1], self.matrix[0][0]))

    def get_translation(self):
        assert(self.matrix.shape == (3, 3))
        return self.matrix[:2,-1:].reshape((2,))
       
def kabsch(crds1, crds2):
    """ returns the rotation matrix to rotate crds2 onto crds1 """
    assert(crds1.shape == crds2.shape)

    correlation_matrix = np.dot(np.transpose(crds1), crds2)
    v, s, w_tr = np.linalg.svd(correlation_matrix)

    is_reflection = (np.linalg.det(v) * np.linalg.det(w_tr)) < 0.0
    if is_reflection:
        s[-1] = -s[-1]
        v[:,-1] = -v[:,-1]
    
    rmat = np.dot(v, w_tr)

    # make sure the matrix is orthonormal 
    det = np.linalg.det(rmat) 
    assert(abs(det - 1.0) < 1e-10)
    
    return rmat

def center_of_mass(coords):
    return coords.mean(axis=0)

class Superposition(object):
    def __init__(self, goal_coords, coords_to_move, transformation):
        self.goal_coords = goal_coords
        self.coords_to_move = coords_to_move
        self.transformation = transformation

    def get_piecewise_rms_before_optimization(self):
        return piecewise_rms(self.goal_coords, self.coords_to_move)

    def get_piecewise_rms(self):
        new_coords = self.transformation.transform(self.coords_to_move)
        return piecewise_rms(self.goal_coords, new_coords)

    def get_euler_angle(self):
        return self.transformation.get_euler_angle()

    def get_translation(self):
        return self.transformation.get_translation()

def superimpose(goal_coords, coords_to_move):
    """ return the transformation to move 'coords_to_move' onto 'goal_coords' """
   
    goal_coords_com = center_of_mass(goal_coords)
    goal_coords_at_origin = goal_coords - goal_coords_com
    
    coords_to_move_com = center_of_mass(coords_to_move)
    coords_to_move_at_origin = coords_to_move - coords_to_move_com

    rmat = kabsch(goal_coords_at_origin, coords_to_move_at_origin)

    transformation = Transformation(translation=-coords_to_move_com)
    transformation.push_front_rotation(rmat)
    transformation.push_front_translation(goal_coords_com)
    return Superposition(goal_coords, coords_to_move, transformation)
