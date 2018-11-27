% 4-10-2018

% turn on the light
ol490.setWhite
    
vid = videoinput('pointgrey', 1, 'F7_Mono8_844x676_Mode5');
src = getselectedsource(vid);

vid.FramesPerTrigger = 1;

preview(vid);

a = input('Press Enter to save location of ROI:')

xy = ludl.getXY
z = ludl.getZ


a = input('Press Enter to save location of white:')

xy_white = ludl.getXY
z_white = ludl.getZ

delete(vid)


