import pygalmesh

p = pygalmesh.Polygon2D([[-0.5, -0.3], [0.5, -0.3], [0.0, 0.5]])
max_edge_size_at_feature_edges = 0.1
domain = pygalmesh.Extrude(
    p,
    [1.0, 0.0, 1.0],
    0.5 * 3.14159265359,
    max_edge_size_at_feature_edges,
)
mesh = pygalmesh.generate_surface_mesh(
    domain,
    min_facet_angle=30.0,
    max_radius_surface_delaunay_ball=0.1,
    max_facet_distance=0.1,
    verbose=False
)

mesh.write("test.obj")
