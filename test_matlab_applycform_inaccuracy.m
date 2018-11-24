% a test program to show that the MATLAB applycform from sRGB to LAB is
% inaccurate
cf = makecform('srgb2lab')
rgb = zeros(1,1,3,'uint8');
rgb(1,1,1:3) = [156 170 56]
rgb(1,1,1:3) = [236 205 0]
lab = applycform(rgb,cf)
labd = lab2double(lab)


cf2 = makecform('lab2srgb')
lab = [82.99057 -0.23356 82.95198]
rgb = applycform(lab,cf2) * 255


