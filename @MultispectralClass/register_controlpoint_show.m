% 11-17-2018
% changed to 6 subplots
% 11-16-2018
% for paper revision
% 8-15-2016
% Register an Aerial Photograph to a Digital Orthophoto
% how to show the aligned images?

function register_controlpoint_show (original_fn, unregistered_fn, folder_output)

%%
% load two images
original = imread(original_fn);
unregistered = imread(unregistered_fn);

% fuse two images with red and green
imf_before = imfuse(original,unregistered,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);

%% Control Points
% use cpselect GUI to choose control points
% and wait until control points are exported

% use previous selected control points?
load([folder_output '/my_control_points'],'movingPoints','fixedPoints');

% cp_method = 'NonreflectiveSimilarity';
cp_method = 'affine';

%% Registered results using control points only
% calculate the transform
mytform_cp = fitgeotrans(movingPoints, fixedPoints, cp_method);
mytform_cp.T

% apply the transform and clip it with the window of the original image
registered_cp = imwarp(unregistered, mytform_cp,'OutputView',imref2d(size(original)));

% fuse two images with red and green
imf_cp = imfuse(original,registered_cp,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);


%% Cross correlation
% refine the control points with cross correlation
% works only for translation (shifting on the XY plane)
% only within 4 pixels: "cpcorr only moves the position of a control point by up to four pixels. Adjusted coordinates are accurate up to one-tenth of a pixel. cpcorr is designed to get subpixel accuracy from the image content and coarse control point selection."
% https://www.mathworks.com/help/images/ref/cpcorr.html?searchHighlight=cpcorr&s_tid=doc_srchtitle

movingPointsAdjusted = cpcorr(movingPoints,fixedPoints,...
    squeeze(unregistered(:,:,2)),...
    squeeze(original(:,:,2)));

% calculate the transform
mytform_cp_corr = fitgeotrans(movingPointsAdjusted, fixedPoints, cp_method);
mytform_cp_corr.T

% apply the transform and clip it with the window of the original image
registered_cp_corr = imwarp(unregistered, mytform_cp_corr,'OutputView',imref2d(size(original)) );

% fuse two images with red and green
imf_cp_corr = imfuse(original,registered_cp_corr,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);

%% save files
save([folder_output '/my_geotrans'],'mytform_cp','mytform_cp_corr')
save([folder_output '/output_images','original'],'unregistered','registered_cp','registered_cp_corr','imf_before','imf_cp_corr')

%% Visualization

figure('Units','inches',...
    'Position',[5 5 6 6],...
    'PaperPositionMode','auto');

set(gca,...
    'Units','normalized',...
    'Position',[.15 .2 .75 .7],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',9,...
    'FontName','Arial');

sc = 100;
sp(1) = subplot(2,2,1);
imshow(original,'InitialMagnification',sc);
title('(a) Truth')
axis image
axis off

sp(2) = subplot(2,2,2);
imshow(unregistered,'InitialMagnification',sc);
title('(b) WSI')
axis image
axis off

sp(3) = subplot(2,2,3);
imshow(imf_before,'InitialMagnification',sc);
title('(c) Images fused before registration')
axis image
axis off

sp(4) = subplot(2,2,4);
imshow(imf_cp_corr,'InitialMagnification',sc);
title('(d) Images fused after registration')
axis image
axis off

if 1
    
    %% tight axis
    % control parameters
    gmarginx = 0.05;
    gmarginy = 0.05;
    gx = 2;
    gy = 2;

    % calculate
    gstepx = (1-gmarginx)/gx;
    gstepy = (1-gmarginy)/gy;
    for i = 1:gx*gy
        grow = (gy-1) - floor((i-1)/gx);
        gcolumn = mod(i-1,gx);
        gposx = gmarginx + gcolumn * gstepx;
        gposy = gmarginy + grow * gstepy;
        sp(i).Position(1) = gposx;
        sp(i).Position(2) = gposy;
        sp(i).Position(3) = gstepx - gmarginx;
        sp(i).Position(4) = (gstepy - gmarginy);
    end
end

return
end
