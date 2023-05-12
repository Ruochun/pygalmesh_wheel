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
def GenWheel(rad=0.25, width=0.2, cp_deviation=0., g_height=0.025, g_width=0.005, g_density=12, g_amp=0., g_period=0, g_curved=False):
    # The angle between 2 adjacent grousers
    g_angle = 2 * math.pi / g_density
    # 1 if this wheel is convex (bump outwards)
    convex = 1 if cp_deviation >= 0. else -1
    # small wiggle room needed in many places
    wiggle_dist = g_width / 2
    small_dist = 0.001

    # the wheel perimeter shape is defined by a 4-point Bezier curve
    outer_CPs = np.empty((num_CPs, 2)).astype(float)
    outer_CPs[0, :] = [rad, width/2]
    outer_CPs[1, :] = [rad + rad*cp_deviation, width/2 - width/3]
    outer_CPs[2, :] = [rad + rad*cp_deviation, width/2 - width/3*2]
    outer_CPs[3, :] = [rad, -width/2]

    inner_CPs = np.array(outer_CPs)
    inner_CPs[:, 0] -= shell_thickness

    # Now create a closed polygon using the 2 perimeters
    num_sample = 20
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
    # rotate it
    wheel_peri = pygalmesh.Rotate(wheel_peri, [1, 0, 0], math.pi / 2)

    # Now add grousers
    # 
    g_seg = []
    g_seed = np.array([0., -width/2, rad + g_height/2 - small_dist]) 
    last_signed_bump = 0.
    # Each grouser period, num of sample points
    num_g_p_sample = 4 if not(g_curved) else 16
    # Start and end t point of each period
    s_e_t = np.linspace(0, 1, g_period+1)
    for i in range(g_period):
        start_t = s_e_t[i]
        end_t = s_e_t[i+1]
        per_seg_lin_length = width / g_period / num_g_p_sample
        period_last_y = 0.
        
        for j in range(1, num_g_p_sample+1):
            t = start_t + (end_t-start_t)/num_g_p_sample*j
            b = pyBernstein(num_CPs - 1, t)
            # eval_pnt[0] gives the vertical move dist
            eval_pnt = b @ outer_CPs
            signed_wheel_bump = eval_pnt[0] - rad
            bump_diff = signed_wheel_bump - last_signed_bump
            last_signed_bump = signed_wheel_bump

            # planal length of this grouser segment
            period_this_y = np.sin(math.pi*2/num_g_p_sample*j)*g_amp
            planal_length = np.sqrt((period_this_y-period_last_y)**2 + per_seg_lin_length**2)
            rot_angle = np.arctan((period_this_y-period_last_y)/per_seg_lin_length)
            period_last_y = period_this_y

            # finally, the cross-section parallelogram shape
            parallelogram_deg = np.arctan(bump_diff/planal_length)
            point2 = [planal_length+np.cos(parallelogram_deg)*wiggle_dist, -bump_diff-np.sin(parallelogram_deg)*wiggle_dist]
            seg_length = np.sqrt(point2[0]**2 + point2[1]**2) - wiggle_dist # smaller... 
            point3 = point2.copy()
            point3[1] += g_height
            parallelogram = pygalmesh.Polygon2D([[0,0], point2, point3, [0, g_height]])

            # create the shape (bugged, only z direction works)
            half1 = pygalmesh.Extrude(
                parallelogram,
                [0.0, 0.0, g_width/2],
                0
            )
            half2 = pygalmesh.Extrude(
                parallelogram,
                [0.0, 0.0, -g_width/2],
                0
            )
            seg = pygalmesh.Union([half1, half2])
            # rotate so facing x
            seg = pygalmesh.Translate(seg, [0, -g_height/2, 0])
            seg = pygalmesh.Rotate(seg, [0,1,0], -math.pi/2)
            seg = pygalmesh.Rotate(seg, [1,0,0], -math.pi/2)
            
            # rotate and move...
            seg = pygalmesh.Rotate(seg, [0,0,1], rot_angle)
            seg = pygalmesh.Translate(seg, g_seed.tolist())
            g_seg.append(seg)
            
            # now where's the new g_seed?
            g_seed += np.array(
                [-seg_length*np.cos(parallelogram_deg)*np.sin(rot_angle),
                seg_length*np.cos(parallelogram_deg)*np.cos(rot_angle),
                seg_length*np.sin(parallelogram_deg)]
                )


    # wheel = pygalmesh.Union(g_seg + [wheel_peri])

    grousers = pygalmesh.Union(g_seg)
    wheel = pygalmesh.Union([grousers, wheel_peri])

    # if geo is disconnected, this function only creates mesh out of one piece of them...
    mesh = pygalmesh.generate_surface_mesh(
        wheel,
        # bounding_sphere_radius = 1.0,
        min_facet_angle = 15.0,
        max_radius_surface_delaunay_ball = 0.001,
        max_facet_distance = 0.001,
        verbose = False
    )

    # Create the grouser

    return mesh

if __name__ == "__main__":
    mesh = GenWheel(g_height=0.025, g_width=0.005, cp_deviation=-0.1, g_density=12, g_amp=0.03, g_period=3, g_curved=False)
    mesh.write("wheel.obj")
