import random
import bpy

def assignMaterials(mesh, k, idx):
    """为每个分割部分分配一个随机的彩色材料"""

    # 清除所有现有的材料
    mesh.materials.clear()

    for i in range(k):
        material = bpy.data.materials.new(''.join(['mat', mesh.name, str(i)]))
        material.diffuse_color = (random.random(), random.random(),
                                  random.random(), 1.0)
        mesh.materials.append(material)

    for i, id in enumerate(idx):
        mesh.polygons[i].material_index = id
