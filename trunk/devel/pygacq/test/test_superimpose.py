from nose.tools import assert_equals

from superimpose import *
from util import rotate
import numpy as np

def test_transformation():
    # test only translation
    trans = np.array([1.0, 2.0])
    translation = Transformation(translation=trans)
    points = np.array([[2.0, 5.0]])

    translated_points = translation.transform(points)
    assert(np.all(translated_points == np.array([[3.0, 7.0]])))

    # test only rotation
    rotation = Transformation(euler_angle=90)
    rotated_points = rotation.transform(points)

    diff = rotated_points - np.array([[-5.0, 2.0]])
    assert(np.all(diff < 1e-10))

    # test translation by (1.0,2.0), then by 90 degree rotation
    both = Transformation()
    both.push_front(translation)
    both.push_front(rotation)
    new_points = both.transform(points)
    diff = new_points - np.array([[-7.0, 3.0]])
    assert(np.all(diff < 1e-10))

    # test 90 degree rotation, then translation by (1.0, 2.0)
    both = Transformation()
    both.push_front(rotation)
    both.push_front(translation)
    new_points = both.transform(points)

    diff = rotated_points - np.array([[-3.0, 7.0]])
    assert(np.all(diff < 1e-10))
    
def test_transformation_rotation():
    for deg in range(-360, 360):
        unit = np.array([1, 0])
        euler_rotation = rotate(unit, deg)

        points = np.array([unit])

        rotation = Transformation(euler_angle=deg)
        matrix_rotation = rotation.transform(points)

        diff = euler_rotation - matrix_rotation
        assert(np.all(diff < 1e-10)), "difference between rotation angle and matrix = %r" % diff

def test_exact_kabsch():
    points1 = np.array([[ 0.0, 0.0],
                        [ 1.0, 0.0],
                        [-1.0, 0.0],
                        [ 0.0, 1.0]])

    for deg in range(-360, 360):
        points2 = np.array(points1)
        rotation = Transformation(euler_angle=deg)
        points2 = rotation.transform(points2)
        rms_before = piecewise_rms(points1, points2)

        # kabsch always returns the shortest rotation
        if deg > 180.0:
            deg = deg - 360
        elif deg < -180:
            deg = deg + 360

        # same direction as rmat
        kabsch_rmat = kabsch(points2, points1)
        kabsch_rotation = Transformation(rotmat=kabsch_rmat)
        rms = piecewise_rms(points2, kabsch_rotation.transform(points1))
        assert(np.all(rms <= rms_before))
        assert(np.all(rms < 1e-10))

        angle = kabsch_rotation.get_euler_angle()
        assert(abs(angle - float(deg)) < 1e-10)

        # opposite direction as rmat
        kabsch_rmat = kabsch(points1, points2)
        kabsch_rotation = Transformation(rotmat=kabsch_rmat)
        rms = piecewise_rms(points1, kabsch_rotation.transform(points2))
        assert(np.all(rms <= rms_before))
        assert(np.all(rms < 1e-10))

        angle = -kabsch_rotation.get_euler_angle()
        assert abs(angle - float(deg)) < 1e-10

def test_nonexact_kabsch():
    """ almost the same as test_exact_kabsch, except now the
    piece-wise rms doesn't necessarily have to be lower if the total
    rmsd is lower."""
    points1 = np.array([[ 0.0, 0.0],
                        [ 1.0, 0.0],
                        [-1.0, 0.0],
                        [ 0.0, 1.0]])

    jittered = np.array([[ 0.1, 0.1],
                         [ 1.1, 0.0],
                         [-1.1, 0.0],
                         [ 0.0, 0.9]])

    for deg in range(-360, 360):
        points2 = np.array(jittered)
        rotation = Transformation(euler_angle=deg)
        points2 = rotation.transform(points2)

        rms_before = piecewise_rms(points1, points2)

        # kabsch always returns the shortest rotation
        if deg > 180.0:
            deg = deg - 360
        elif deg < -180:
            deg = deg + 360

        # same direction as rmat
        kabsch_rmat = kabsch(points2, points1)
        kabsch_rotation = Transformation(rotmat=kabsch_rmat)
        rms = piecewise_rms(points2, kabsch_rotation.transform(points1))
        assert(rmsd(rms) <= rmsd(rms_before))

        angle = kabsch_rotation.get_euler_angle()
        assert abs(angle - float(deg)) < 1e-10

        # opposite direction as rmat
        kabsch_rmat = kabsch(points1, points2)
        kabsch_rotation = Transformation(rotmat=kabsch_rmat)
        rms = piecewise_rms(points1, kabsch_rotation.transform(points2))
        assert(rmsd(rms) <= rmsd(rms_before))

        angle = -kabsch_rotation.get_euler_angle()
        assert abs(angle - float(deg)) < 1e-10

test_nonexact_kabsch()

def test_kabsch_with_translation():
    """ note, the angles are no longer reproducible, yet the RMSD should still be smaller """
    points1 = np.array([[ 0.0, 0.0],
                        [ 1.0, 0.0],
                        [-1.0, 0.0],
                        [ 0.0, 1.0]])

    for deg in range(-360, 360):
        points2 = points1 + np.array([1.0, 0.0])

        rotation = Transformation(euler_angle=deg)
        points2 = rotation.transform(points2)
        rms_before = piecewise_rms(points1, points2)

        # kabsch always returns the shortest rotation
        if deg > 180.0:
            deg = deg - 360
        elif deg < -180:
            deg = deg + 360

        # same direction as rmat
        kabsch_rmat = kabsch(points2, points1)
        kabsch_rotation = Transformation(rotmat=kabsch_rmat)
        rms = piecewise_rms(points2, kabsch_rotation.transform(points1))
        assert(rmsd(rms) <= rmsd(rms_before))

        # opposite direction as rmat
        kabsch_rmat = kabsch(points1, points2)
        kabsch_rotation = Transformation(rotmat=kabsch_rmat)
        rms = piecewise_rms(points1, kabsch_rotation.transform(points2))
        assert(rmsd(rms) <= rmsd(rms_before))

def assert_superimpose_successful(points1, points2, expected_rmsd=None):
    superposition = superimpose(points1, points2)

    rms_before = superposition.get_piecewise_rms_before_optimization()
    rms_after = superposition.get_piecewise_rms()

    assert(rmsd(rms_after) <= rmsd(rms_before))
    if expected_rmsd is not None:
        assert(rmsd(rms_after) <= expected_rmsd)

    return superposition

def test_superimpose():
    points1 = np.array([[ 0.0, 0.0],
                        [ 1.0, 0.0],
                        [-1.0, 0.0],
                        [ 0.0, 1.0]])

    points2 = points1 + np.array([1.0, 0.0])
    assert_superimpose_successful(points1, points2, 0.0)

    points1 += np.array([1.2345, 6.789])
    assert_superimpose_successful(points1, points2, 1e-10)

    points2 -= np.array([9.8765, 4.321])
    assert_superimpose_successful(points1, points2, 1e-10)

    for trans in range(10):
        for deg in range(-360, 360):
            rpoints = points2 + np.array([float(trans), float(trans)])
            rotation = Transformation(euler_angle=deg)
            rpoints = rotation.transform(rpoints)
            assert_superimpose_successful(points1, rpoints, 1e-10)

    # some tests with jittered points
    jittered = np.array([[ 0.1,-0.1],
                         [ 1.1, 0.0],
                         [-1.1, 0.0],
                         [ 0.0, 0.9]])
    assert_superimpose_successful(points1, jittered, 0.12)

    jittered += np.array([3.1415, 9.265])
    assert_superimpose_successful(points1, jittered, 0.12)

    for trans in range(10):
        for deg in range(-360, 360):
            rpoints = jittered + np.array([float(trans), float(trans)])
            rotation = Transformation(euler_angle=deg)
            rpoints = rotation.transform(rpoints)
            assert_superimpose_successful(points1, rpoints, 0.12)


def test_superimpose_with_real_data():

    box_centers = np.array([(1330.0, 4220.0),
                            (2708.5, 3518.5),
                            (2990.0,  642.0)])
    
    star_centers = np.array([(1325.69, 4221.12),
                             (2704.22, 3518.40),
                             (2984.30,  641.47)])

    field_center = np.array([3109, 2304])
    box_centers -= field_center
    star_centers -= field_center

    superposition = assert_superimpose_successful(star_centers, box_centers, 0.48)

    rms_before = superposition.get_piecewise_rms_before_optimization()
    rms_after  = superposition.get_piecewise_rms()
    print rmsd(rms_before), "->", rmsd(rms_after)
    print rms_before
    print rms_after
    print superposition.get_euler_angle()
    print superposition.get_translation()
