import vtk

# Read from file
stlreader = vtk.vtkSTLReader()
stlreader.SetFileName("dental_0173.stl")

# 输出所有点坐标
stlreader.Update()
stlMapper = vtk.vtkPolyDataMapper()
stlMapper.SetInputConnection(stlreader.GetOutputPort())
p = [0, 0, 0]
polydata = stlreader.GetOutput()
for i in range(polydata.GetNumberOfPoints()):
    polydata.GetPoint(i, p)
    print(p)

cylinderMapper = vtk.vtkPolyDataMapper()  # 渲染多边形几何数据
cylinderMapper.SetInputConnection(stlreader.GetOutputPort())  # VTK可视化管线的输入数据接口 ，对应的可视化管线输出数据的接口为GetOutputPort()；
cylinderActor = vtk.vtkActor()
cylinderActor.SetMapper(cylinderMapper)  # 设置生成几何图元的Mapper。即连接一个Actor到可视化管线的末端(可视化管线的末端就是Mapper)。
renderer = vtk.vtkRenderer()  # 负责管理场景的渲染过程
renderer.AddActor(cylinderActor)
renderer.SetBackground(0.1, 0.2, 0.4)
renWin = vtk.vtkRenderWindow()  # 将操作系统与VTK渲染引擎连接到一起。
renWin.AddRenderer(renderer)
renWin.SetSize(300, 300)
iren = vtk.vtkRenderWindowInteractor()  # 提供平台独立的响应鼠标、键盘和时钟事件的交互机制
iren.SetRenderWindow(renWin)

# 交互器样式的一种，该样式下，用户是通过控制相机对物体作旋转、放大、缩小等操作
style = vtk.vtkInteractorStyleTrackballCamera()

iren.SetInteractorStyle(style)
iren.Initialize()
iren.Start()

# Clean up
# del cylinder
del stlreader
del cylinderMapper
del cylinderActor
del renderer
del renWin
del iren