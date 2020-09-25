function XYZ = Yxy2XYZ (Yxy)
    Y = Yxy(1);
    x = Yxy(2);
    y = Yxy(3);
    z = 1 - x - y;
    X = x*Y/y;
    Z = z*Y/y;
    XYZ = [X Y Z];
end

