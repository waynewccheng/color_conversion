Yxy = [501.81 0.3995 0.3975];
XYZ = Yxy2XYZ(Yxy)

Yxy0 = [1163.82 0.3941 0.3892];
XYZ0 = Yxy2XYZ(Yxy0)

lab1 = XYZ2lab(XYZ,XYZ0)

Yxy = [56.73 0.2922 0.3123];
XYZ = Yxy2XYZ(Yxy)

Yxy0 = [128.41 0.2888 0.3072];
XYZ0 = Yxy2XYZ(Yxy0)

lab2 = XYZ2lab(XYZ,XYZ0)

dE = lab2de00(lab1,lab2)
