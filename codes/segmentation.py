import bpy
import mathutils
import math
import numpy
import scipy.linalg
import scipy.cluster
import scipy.sparse
import scipy.sparse.csgraph
import scipy.sparse.linalg


delta = None


eta = None


def _face_center(mesh, face):
    """计算给定面中心的坐标"""
    center = mathutils.Vector()
    for vert in face.vertices:
        center += mesh.vertices[vert].co
    return center/len(face.vertices)


def _geodesic_distance(mesh, face1, face2, edge):
    """计算两个相邻面face1和face2之间给定边上的测地线距离"""
    edge_center = (mesh.vertices[edge[0]].co + mesh.vertices[edge[1]].co)/2
    return (edge_center - _face_center(mesh, face1)).length + \
           (edge_center - _face_center(mesh, face2)).length


def _angular_distance(mesh, face1, face2):
    """计算给定相邻面的角距离"""
    angular_distance = (1 - math.cos(mathutils.Vector.angle(face1.normal,
                                                            face2.normal)))
    if (face1.normal.dot(_face_center(mesh, face2) -
                         _face_center(mesh, face1))) < 0:

        angular_distance *= eta
    return angular_distance


def _create_distance_matrix(mesh):
    """创建所有相邻面之间的角距离和测地线距离矩阵。返回矩阵的第i,j项包含第i个面和第j个面之间的距离."""

    faces = mesh.polygons
    l = len(faces)

    # 从边键映射到相邻的面
    adj_faces_map = {}
    # 通过迭代边来找到相邻的面
    for index, face in enumerate(faces):
        for edge in face.edge_keys:
            if edge in adj_faces_map:
                adj_faces_map[edge].append(index)
            else:
                adj_faces_map[edge] = [index]

    # 帮助向量创建稀疏矩阵
    row_indices = []
    col_indices = []
    Gval = []  # 角距离矩阵的值
    Aval = []  # 测地线距离矩阵的值
    # 迭代相邻面并计算距离
    for edge, adj_faces in adj_faces_map.items():
        if len(adj_faces) == 2:
            i = adj_faces[0]
            j = adj_faces[1]

            Gtemp = _geodesic_distance(mesh, faces[i], faces[j], edge)
            Atemp = _angular_distance(mesh, faces[i], faces[j])
            Gval.append(Gtemp)
            Aval.append(Atemp)
            row_indices.append(i)
            col_indices.append(j)
            
            Gval.append(Gtemp)
            Aval.append(Atemp)
            row_indices.append(j)
            col_indices.append(i)

        elif len(adj_faces) > 2:
            print("Edge with more than 2 adjacent faces: " + str(adj_faces) + "!")

    Gval = numpy.array(Gval)
    Aval = numpy.array(Aval)
    values = delta * Gval / numpy.mean(Gval) + \
             (1.0 - delta) * Aval / numpy.mean(Aval)

    # 创建稀疏矩阵
    distance_matrix = scipy.sparse.csr_matrix(
        (values, (row_indices, col_indices)), shape=(l, l))
    return distance_matrix


def _create_affinity_matrix(mesh):
    """创建给定网格的邻接矩阵"""

    l = len(mesh.polygons)
    print("mesh_segmentation: Creating distance matrices...")
    distance_matrix = _create_distance_matrix(mesh)

    print("mesh_segmentation: Finding shortest paths between all faces...")
    
    W = scipy.sparse.csgraph.dijkstra(distance_matrix)
    inf_indices = numpy.where(numpy.isinf(W))
    W[inf_indices] = 0

    print("mesh_segmentation: Creating affinity matrix...")
    # 将距离更改为相似点
    sigma = W.sum()/(l ** 2)
    den = 2 * (sigma ** 2)
    W = numpy.exp(-W/den)
    W[inf_indices] = 0
    numpy.fill_diagonal(W, 1)

    return W


def _initial_guess(Q, k):
    """计算簇中心的初始猜测，以贪婪的方式选择彼此关联最小的观测指标。Q是观测值的关联矩阵."""

    # 选择彼此关联程度最低的一对索引
    min_indices = numpy.unravel_index(numpy.argmin(Q), Q.shape)

    chosen = [min_indices[0], min_indices[1]]
    for _ in range(2,k):

        new_index = numpy.argmin(numpy.max(Q[chosen,:], axis=0))
        chosen.append(new_index)

    return chosen


def segment_mesh(mesh, k, coefficients, action, ev_method, kmeans_init):
    """将给定的网格划分为k个簇，并对每个簇执行给定的操作"""

    # 设置参数
    global delta
    global eta
    delta, eta = coefficients

    # 关联矩阵
    W = _create_affinity_matrix(mesh)
    print("mesh_segmentation: Calculating graph laplacian...")
    # 次数矩阵
    Dsqrt = numpy.sqrt(numpy.reciprocal(W.sum(1)))
    # 图拉普拉斯算子
    L = ((W * Dsqrt).transpose() * Dsqrt).transpose()

    print("mesh_segmentation: Calculating eigenvectors...")
    # 得到特征向量
    if ev_method == 'dense':
        _, V = scipy.linalg.eigh(L, eigvals = (L.shape[0] - k, L.shape[0] - 1))
    else:
        _, V = scipy.sparse.linalg.eigsh(L, k)
    # 将每一行规范化为单位长度
    V /= numpy.linalg.norm(V, axis=1)[:,None]

    if kmeans_init == 'kmeans++':
        print("mesh_segmentation: Applying kmeans...")
        _, idx = scipy.cluster.vq.kmeans2(V, k, minit='++', iter=50)
    else:
        print("mesh_segmentation: Preparing kmeans...")
        # 计算关联矩阵
        Q = V.dot(V.transpose())
        # 计算聚类的初始值
        initial_centroids = _initial_guess(Q, k)

        print("mesh_segmentation: Applying kmeans...")
        _, idx = scipy.cluster.vq.kmeans2(V, V[initial_centroids,:], iter=50)

    print("mesh_segmentation: Done clustering!")
    # 对集群结果执行操作
    if action:
        action(mesh, k, idx)
