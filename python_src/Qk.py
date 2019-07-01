import phom

#voro: phom.Voronoi2D object

#vector from particle centers to centroids, flattened (cx1, cy1, cx2, cy2, ....)
def c_vec(voro):
    centroids = voro.get_cell_centroids()
    centers = np.dot((voro.embed.box_mat, voro.embed.pos.reshape(2,-1))).transpose().flatten()

    return centroids - centers


#given a vector field defined on the particles and the indices of a triangle, compute discrete divergence
def compute_triangle_divergence(field, indices, voro):
    full_indices = np.concatenate([[2*i, 2*i+1] for i in indices])
    RHS = field[full_indices]

    pos = np.dot((voro.embed.box_mat, voro.embed.pos[full_indices].reshape(2,-1))).transpose()
    mat = []
    mat.append([1, 0, pos[0,0], pos[0,1], 0, 0])
    mat.append([1, 0, 0, 0 pos[0,0], pos[0,1]])
    mat.append([1, 0, pos[1,0], pos[1,1], 0, 0])
    mat.append([1, 0, 0, 0 pos[1,0], pos[1,1]])
    mat.append([1, 0, pos[2,0], pos[2,1], 0, 0])
    mat.append([1, 0, 0, 0 pos[2,0], pos[2,1]])

    mat = np.concatenate(mat)

    d = np.linalg.solve(mat, RHS)

    return d[2] + d[5]
