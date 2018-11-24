function lab = srgb2lab (rgb)
    xyz = srgb2xyz(rgb);
    rgbn = [1 1 1];
    xyzn = srgb2xyz(rgbn);

    lab = xyz2lab(xyz(1),xyz(2),xyz(3),  xyzn(1),xyzn(2),xyzn(3));
end

function rgb = xyY2srgb (xyY)
% http://en.wikipedia.org/wiki/SRGB
    mat = [3.2306 -1.5372 -0.4986; -0.9689 1.8758 0.0415; 0.0557 -0.2040 1.0570];

    % calculate tristimulus
    x = xyY(1);
    y = xyY(2);
    Y = xyY(3);
    z = 1 - x - y;
    X = x * (Y/y);
    Z = z * (Y/y);
    
    % transform
    rgblin = mat * [X Y Z]';
    rgblin = min(rgblin,1);
    rgblin = max(rgblin,0);

    % linearize
    rgb = srgblin(rgblin);

    return
        
    % helper function for xyz2srgb
    function c_rgb = srgblin (clin)
        a = 0.055;
        if clin <= 0.0031308
            c_rgb = clin .* 12.92;
        else
            c_rgb = power(clin,1/2.4).*(1+a) - a; 
        end
        c_rgb = min(c_rgb,1);
        c_rgb = max(c_rgb,0);
    end
end

function ret = xyz2lab (x, y, z, xn, yn, zn)
retl = 116 * f3(y/yn) - 16;
reta = 500 * (f3(x/xn) - f3(y/yn));
retb = 200 * (f3(y/yn) - f3(z/zn));
ret = [retl reta retb];
end

function ret = f3 (t)
if t > power(6/29,3)
    ret = power(t,1/3);
else
    ret = t/3*power(29/6,2)+4/29;
end
end


function XYZxyz = srgb2xyz (rgb)
    mat = [0.4124 0.3576 0.1805; 0.2126 0.7152 0.0722; 0.0193 0.1192 0.9505];
    rgblin = srgblin(rgb);
    XYZ = mat * rgblin';
    sumXYZ = XYZ(1) + XYZ(2) + XYZ(3);
    if (sumXYZ == 0)
        xyz = [1/3 1/3 1/3];
    else
        xyz = XYZ ./ sumXYZ;
    end
    XYZxyz = [XYZ' xyz'];
    
    % helper function for srgv2xyz
    function c_linear = srgblin (csrgb)
        a = 0.055;
        if csrgb <= 0.04045
            c_linear = csrgb ./ 12.92;
        else
            c_linear = power((csrgb+a)./(1+a),2.4); 
        end
    end
end
