% 10-22-2013
% 
% 7-31-2013
% fixed RGBScan plot problem by adding sf1, sf2, and sf3
% input: image file, display data
% output: imlab1=original imlab2=display [n*3]

function image2lab (imagename, displayname)
    global mea;
    global specarray;

        
    % load display data
    % displayname = 'BARCO_GMA22_WOVD_sRGB_10-05-2012_1638';
    
datafile = ['display/' displayname '.mat']

%    datafile = 'display/cinema_11-23-2010_1633.mat';
%datafile = 'dell19r_11-23-2010_1852.mat';
%datafile = 'chimei_11-23-2010_2022.mat';
%datafile = 'iphone4_11-24-2010_1252.mat';
%datafile = 'nexus1_12-20-2010_1525.mat';
%datafile = '3dvision_11-23-2010_1102.mat';
%datafile = 'galaxys_12-20-2010_1252.mat';
%datafile = 'r31_11-23-2010_1352.mat';
%datafile = 'nec_pa271w_03-01-2011_1625.mat';
%datafile = 'display/SONY PVM_2551MD_60Hz_09-08-2012_0625'
%datafile = 'display/DELL 2212_08-21-2012_1155'
%datafile = 'r31_07-31-2013_0759'
%datafile = 'Dell19_07-31-2013_0614'
%datafile = 'Samsung19_08-01-2013_1156'
%datafile = 'NEC_PA271W_08-02-2013_0701'


load(datafile);

    % white point
    % for reference white in CIELAB
    wY = mea(1024,4);
    wx = mea(1024,5);
    wy = mea(1024,6);
    wz = 1 - wx - wy;
    wX = wx * (wY/wy);
    wZ = wz * (wY/wy);
    [wX wY wZ wx wy wz]
    
    % black point
    % for display model: R+G+B-2K   
    kY = mea(1,4);
    kx = mea(1,5);
    ky = mea(1,6);
    kz = 1 - kx - ky;
    kX = kx * (kY/ky);
    kZ = kz * (kY/ky);


    % load image
    im = imread(imagename);
    [dimy dimx dimz] = size(im);


    % calculate input image
    imlabinput = zeros(dimy,dimx,3);
    for j = 1:dimy
        for i = 1:dimx
            rgb = squeeze(double(im(j,i,1:3)))/255;
            [lstar astar bstar] = srgb2lab(rgb');
            imlabinput(j,i,1:3) = [lstar astar bstar] ;
        end
    end
    imlab1 = reshape(imlabinput,dimx*dimy,3);
    
    
    % calculate output image
    % causion type casting; spent 6 hours in debugging
    imlaboutput = zeros(dimy,dimx,3);
    imrgboutput = zeros(dimy,dimx,3);
    
    for j = 1:dimy
        for i = 1:dimx
            r = double(im(j,i,1));
            g = double(im(j,i,2));
            b = double(im(j,i,3));
            [lstar astar bstar] = xyzpredict(r, g, b);
            imlaboutput(j,i,1:3) = [lstar astar bstar] ;
            imrgboutput(j,i,1:3) = lab2srgb([lstar astar bstar]);
        end
    end
    imlab2 = reshape(imlaboutput,dimx*dimy,3);

    save('imlab','imlab1','imlab2')

    imrgboutput8 = im;
    imrgboutput8 = uint8(imrgboutput);
    
    imagenamehead = imagename(1:strfind(imagename,'.')-1);
    imwrite(imrgboutput8,[imagenamehead '_' displayname '.jpg'],'jpg')

    
    %---------------------------- visualization
    clf

    % show original image
    subplot(3,3,1)
    image(im)
    axis ij
    axis off
    title('Image')
    
    % show L a b as curves
     subplot(3,3,2)
     image(imrgboutput/255)
     
     
%     hold on
%     plot(imlab2(:,1),'r')
%     plot(imlab2(:,2),'g')
%     plot(imlab2(:,3),'b')
% 


    % scale factor; by tril
    sf1 = 100;
    sf2 = 300;
    sf3 = 300;


    % show before L a b
    subplot(3,3,4)
    k1 = double(imlaboutput(:,:,1)/sf1);
    image(repmat(k1,[1,1,3]))
    axis ij
    axis off
    title('Input L*')    

    subplot(3,3,5)
    k2 = double(imlaboutput(:,:,2)/sf2 + 0.5);
    image(repmat(k2,[1,1,3]))
    axis ij
    axis off
    title('Input a*')
    
    subplot(3,3,6)
    k3 = double(imlaboutput(:,:,3)/sf3 + 0.5);
    image(repmat(k3,[1,1,3]))
    axis ij
    axis off
    title('Input b*')

    % show after L a b
    
    subplot(3,3,7)
    kk1 = double(imlabinput(:,:,1)/sf1);
    image(repmat(kk1,[1,1,3]))
    axis ij
    axis off
    title('Output L*')
    
    subplot(3,3,8)
    kk2 = double(imlabinput(:,:,2)/sf2 + 0.5);
    image(repmat(kk2,[1,1,3]))
    axis ij
    axis off
    title('Output a*')
    
    subplot(3,3,9)
    kk3 = double(imlabinput(:,:,3)/sf3 + 0.5);
    image(repmat(kk3,[1,1,3]))
    axis ij
    axis off
    title('Output b*')
    
    % show difference map
    subplot(3,3,3)
    kk3 = double(imlabinput(:,:,3)/sf3 + 0.5);
    cm3 = abs(kk3 - k3);
    %image(repmat(kk3-k3,[1,1,3]))
    mesh(cm3)
    axis ij

    
    return

    
    % ----------------------- helper functions
    % a display model based on RGBW characterization data
    % input: DDL r, g, b
    % input: XYZ of black point; kX kY kZ from parent function
    % input: XYZ of white point; wX wY wZ from parent function
    % output: L* a* b*
    % method: r+g+b-2k 
    function [lstar astar bstar] = xyzpredict (r, g, b)
        % retrieve XYZ of red channel
        rY = mea(r*4+1,4);
        rx = mea(r*4+1,5);
        ry = mea(r*4+1,6);
        rz = 1 - rx - ry;
        rX = rx * (rY/ry);
        rZ = rz * (rY/ry);

        % retrieve XYZ of green channel
        gY = mea(g*4+2,4);
        gx = mea(g*4+2,5);
        gy = mea(g*4+2,6);
        gz = 1 - gx - gy;
        gX = gx * (gY/gy);
        gZ = gz * (gY/gy);
        
        % retrieve XYZ of blue channel
        bY = mea(b*4+3,4);
        bx = mea(b*4+3,5);
        by = mea(b*4+3,6);
        bz = 1 - bx - by;
        bX = bx * (bY/by);
        bZ = bz * (bY/by);

        % plug in model
        % and convert into L* a* b*
        [lstar astar bstar] = xyz2lab(rX+gX+bX-2*kX, rY+gY+bY-2*kY, rZ+gZ+bZ-2*kZ, wX*1.0, wY*1.0, wZ*1.0);
    end
    
end

function XYZ = srgbpredict (rgb)
        rY = 0.2126;
        rx = 0.6400;
        ry = 0.3300;
        rz = 1 - rx - ry;
        rX = rx * (rY/ry);
        rZ = rz * (rY/ry);
        
        gY = 0.7153;
        gx = 0.3000;
        gy = 0.6000;
        gz = 1 - gx - gy;
        gX = gx * (gY/gy);
        gZ = gz * (gY/gy);
        
        bY = 0.0721;
        bx = 0.1500;
        by = 0.0600;
        bz = 1 - bx - by;
        bX = bx * (bY/by);
        bZ = bz * (bY/by);
    
        mat = [rX gX bX; rY gY bY; rZ gZ bZ];    
        rgblin = srgblin(rgb);
        XYZ = mat * rgblin';
        
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

function [lstar astar bstar] = srgb2lab (rgb)
    XYZxyz = srgb2xyz(rgb);
    XYZxyzn = srgb2xyz([1 1 1]);
    [lstar astar bstar] = xyz2lab (XYZxyz(1), XYZxyz(2), XYZxyz(3), XYZxyzn(1), XYZxyzn(2), XYZxyzn(3));
end

function [lstar astar bstar] = xyz2lab (x, y, z, xn, yn, zn)
    lstar = 116 * helpf(y/yn) - 16;
    astar = 500 * (helpf(x/xn) - helpf(y/yn));
    bstar = 200 * (helpf(y/yn) - helpf(z/zn));

    if lstar > 100
        [x y z xn yn zn lstar astar bstar]
        lstar = 100;
    end
    
    return
    
    function ys = helpf (t)
        if t > power(6/29,3)
            ys = power(t,1/3);
        else
            ys = power(29/6,2)/3*t+4/29;
        end
    end
end

function XYZxyz = srgb2xyz (rgb)
    mat = [0.4124 0.3576 0.1805; 0.2126 0.7152 0.0722; 0.0193 0.1192 0.9505];
    rgblin = srgblin(rgb);
    XYZ = mat * rgblin';
    sumXYZ = XYZ(1) + XYZ(2) + XYZ(3);
    if (sumXYZ == 0)
        xyz = [1/3 1/3 1/3]';
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

function XYZ = xyY2XYZ (x, y, Y)
        z = 1 - x - y;
        X = x * (Y/y);
        Z = z * (Y/y);
        XYZ = [X Y Z];
end
