import pygalmesh
import math
import numpy as np
from scipy.special import comb

shell_thickness = 0.02
max_edge_size_at_feature_edges = 0.1
num_CPs = 4

def pyBernstein(degree, t):
    B = np.empty(degree + 1).astype(float)
    for i in range(len(B)):
        B[i] = comb(degree, i)*(t**i)*(1.0 - t)**(degree - i)
    return B

# 2D shapes are in x-z plane. After 3D wheels are created, they are rotated about x to make them face x-forward.
# cp_deviation is the percentage that the mid 2 control points defining the wheel perimeter deviates from the position where the wheel surface is perfect flat, in z direction 
def GenWheel(rad=0.25, width=0.2, cp_deviation=0., mu=0.4, g_height=0.025, g_width=0.005, g_density=12, g_amp=0., g_period=0, g_curved=False):
    # The angle between 2 adjacent grousers
    g_angle = 2 * math.pi / g_density

    # the wheel perimeter shape is defined by a 4-point Bezier curve
    outer_CPs = np.empty((num_CPs, 2)).astype(float)
    outer_CPs[0, :] = [rad, width/2]
    outer_CPs[1, :] = [rad + rad*cp_deviation, width/2 - width/3]
    outer_CPs[2, :] = [rad + rad*cp_deviation, width/2 - width/3*2]
    outer_CPs[3, :] = [rad, -width/2]

    inner_CPs = np.array(outer_CPs)
    inner_CPs[:, 0] -= shell_thickness

    # Now create a closed polygon using the 2 perimeters
    num_sample = 10
    bz_t = np.linspace(0, 1, num_sample)
    B = np.zeros((num_sample, num_CPs)).astype(float)
    for i in range(num_sample):
        # degree 3 beziers
        B[i, :] = pyBernstein(num_CPs - 1, bz_t[i])
    # created sampled point list for outer and inner perimeters
    outer_list = B @ outer_CPs
    inner_list = B @ inner_CPs
    # Ask pygalmesh to create the shape
    wheel_cross_sec = pygalmesh.Polygon2D(outer_list.tolist() + (np.flip(inner_list, axis=0)).tolist())
    wheel_peri = pygalmesh.RingExtrude(wheel_cross_sec, max_edge_size_at_feature_edges)



    mesh = pygalmesh.generate_surface_mesh(
        wheel_peri,
        # bounding_sphere_radius = 1.0,
        min_facet_angle = 15.0,
        max_radius_surface_delaunay_ball = 0.002,
        max_facet_distance = 0.002,
        verbose = False
    )

    # Create the grouser

    return mesh

if __name__ == "__main__":
    mesh = GenWheel()
    mesh.write("wheel.obj")
