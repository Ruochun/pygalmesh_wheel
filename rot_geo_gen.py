import pygalmesh

p = pygalmesh.Polygon2D([[0.5, -0.3], [1.5, -0.3], [1.0, 0.5]])
max_edge_size_at_feature_edges = 0.1
domain = pygalmesh.RingExtrude(p, max_edge_size_at_feature_edges)
mesh = pygalmesh.generate_surface_mesh(
    domain,
    bounding_sphere_radius = 0.0,
    min_facet_angle = 3.0,
    max_radius_surface_delaunay_ball = 0.1,
    max_facet_distance = 0.1,
    verbose = False,
)

def GenWheel(rad=0.25, width=0.2, curvature=0., mu=0.4, g_height=0.025, g_density=12, g_amp=0., g_period=0, g_curved=False, filename):

    return None

mesh.write("out.obj")
