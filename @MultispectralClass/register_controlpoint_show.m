% 11-17-2018
% changed to 6 subplots
% 11-16-2018
% for paper revision
% 8-15-2016
% Register an Aerial Photograph to a Digital Orthophoto
% how to show the aligned images?

function register_controlpoint_show (original_fn, unregistered_fn, folder_output, lab_truth_fn)

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
if 0
    
    load([folder_output '/my_control_points'],'movingPoints','fixedPoints');
    
    % cp_method = 'NonreflectiveSimilarity';
    cp_method = 'affine';
    
    %% Registered results using control points only
    % calculate the transform
    mytform_cp = fitgeotrans(movingPoints, fixedPoints, cp_method);
    mytform_cp.T
    
    % apply the transform and clip it with the window of the original image
    registered_cp = imwarp(original, mytform_cp,'OutputView',imref2d(size(unregistered)));
    
    % fuse two images with red and green
    imf_cp = imfuse(registered_cp,unregistered,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
   
    %% Cross correlation
    % refine the control points with cross correlation
    % works only for translation (shifting on the XY plane)
    % only within 4 pixels: "cpcorr only moves the position of a control point by up to four pixels. Adjusted coordinates are accurate up to one-tenth of a pixel. cpcorr is designed to get subpixel accuracy from the image content and coarse control point selection."
    % https://www.mathworks.com/help/images/ref/cpcorr.html?searchHighlight=cpcorr&s_tid=doc_srchtitle
    
    movingPointsAdjusted = cpcorr(movingPoints,fixedPoints,...
        squeeze(original(:,:,2)),...
        squeeze(unregistered(:,:,2)));
    
    % calculate the transform
    mytform_cp_corr = fitgeotrans(movingPointsAdjusted, fixedPoints, cp_method);
    mytform_cp_corr.T
    
    % apply the transform and clip it with the window of the original image
    registered_cp_corr = imwarp(original, mytform_cp_corr,'OutputView',imref2d(size(unregistered)) );
    
    % fuse two images with red and green
    imf_cp_corr = imfuse(registered_cp_corr,unregistered,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    
    %% save files
    whos
    save([folder_output '/my_geotrans'],'mytform_cp','mytform_cp_corr')
    save([folder_output '/output_images'],'original','unregistered','registered_cp','registered_cp_corr','imf_before','imf_cp_corr')
    
else
    
    % use Matlab app
    % https://www.mathworks.com/help/images/register-images-using-the-registration-estimator-app.html
    load([folder_output '/matlab_registration'],'movingReg');

    mytform_cp_corr = movingReg.Transformation;
    
    % apply the transform and clip it with the window of the original image
    registered_cp_corr = imwarp(original, mytform_cp_corr,'OutputView',imref2d(size(unregistered)) );
    
    % fuse two images with red and green
    imf_cp_corr = imfuse(registered_cp_corr,unregistered,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    
end


%% trim images
truth_sizex = size(original,2);
truth_sizey = size(original,1);

p_11 = [1 1 1] * mytform_cp_corr.T;
p_x1 = [truth_sizex 1 1] * mytform_cp_corr.T;
p_1y = [1 truth_sizey 1] * mytform_cp_corr.T;
p_xy = [truth_sizex truth_sizey 1] * mytform_cp_corr.T;

P_x_min = ceil(max([1 p_11(1) p_1y(1)]));
P_x_max = floor(min([size(registered_cp_corr,2) p_x1(1) p_xy(1)]));

P_y_min = ceil(max([1 p_11(2) p_x1(2)]));
P_y_max = floor(min([size(registered_cp_corr,1) p_1y(2) p_xy(2)]));

[P_x_min P_x_max P_y_min P_y_max]

unregistered_trimmed = unregistered(P_y_min:P_y_max,P_x_min:P_x_max,:);
registered_cp_corr_trimmed = registered_cp_corr(P_y_min:P_y_max,P_x_min:P_x_max,:);


%% LAB of truth

% load the truth LAB data 
load([lab_truth_fn '/LAB.mat'],'LAB')

% use the same transform to process the LAB data
LAB_truth_reg = imwarp(LAB, mytform_cp_corr,'OutputView',imref2d(size(unregistered)) );

% trim the LAB data
LAB_truth_reg_trimmed = LAB_truth_reg(P_y_min:P_y_max,P_x_min:P_x_max,:);


%% LAB of scan
LAB_scan = rgb2lab(unregistered_trimmed,'ColorSpace','srgb','WhitePoint','d65');


%% save the LAB data
save([folder_output '/final_images'],'unregistered_trimmed','registered_cp_corr_trimmed','LAB_truth_reg_trimmed','LAB_scan')


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
