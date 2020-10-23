% ACL 10/9/20
% Function to convert data output by Paul's HIMS code to a new format
%
% This script assumes that you have Transmittance, XYZ, LAB, & sRGB data from the HIMS scanner
% It will create the necessary folders & convert the necessary files to run
% the MultispectralClass code to register images and compute dE
%
% It will require inputs of:
% path -- the general path in which your data will be saved
% path_processed -- the general path for the data that is to be converted
% Sample_Name -- the name of the sample of interest
% Scanner 1 -- the name of the scanner that is being compared (lowercase)
% Scanner 2 -- the name of the other scanner that is being compared(lowercase)
% Date -- the date the HIMS data was collected (in the form MMDDYY)

function fileconverter(path, path_processed, Sample_Name, Scanner1, Scanner2, Scanner3, Date)

% Make truth directories for MultispectralClass 
mkdir([path '/' Sample_Name '/truth/220 transmittance']);
mkdir([path '/' Sample_Name '/truth/420 CIEXYZ']);
mkdir([path '/' Sample_Name '/truth/520 CIELAB']);
mkdir([path '/' Sample_Name '/truth/900 sRGB']);

% Make scanner directories for Multispectral Class
mkdir([path '/' Sample_Name '/' Scanner1 '/100 wsi']);
mkdir([path '/' Sample_Name '/' Scanner1 '/200 roi']);
mkdir([path '/' Sample_Name '/' Scanner1 '/400 sRGB']);

mkdir([path '/' Sample_Name '/' Scanner2 '/100 wsi']);
mkdir([path '/' Sample_Name '/' Scanner2 '/200 roi']);
mkdir([path '/' Sample_Name '/' Scanner2 '/400 sRGB']);

mkdir([path '/' Sample_Name '/' Scanner3 '/100 wsi']);
mkdir([path '/' Sample_Name '/' Scanner3 '/200 roi']);
mkdir([path '/' Sample_Name '/' Scanner3 '/400 sRGB']);

% Load Transmittance, XYZ, LAB, RGB, and truth image data output by Paul's HIMS code
load([path_processed '/' Date '/' Sample_Name '/Transmittance/trans_mean_camera']);
load([path_processed '/' Date '/' Sample_Name '/CIE_Coord/XYZ_array.mat'],'XYZ_array');
load([path_processed '/' Date '/' Sample_Name '/CIE_Coord/LAB_array.mat'],'LAB_array');
load([path_processed '/' Date '/' Sample_Name '/EndResults/rgb.mat']);
truth = imread([path_processed '/' Date '/' Sample_Name '/EndResults/truth.tif']);

% Reshape LAB and XYZ into 3D arrays
LAB = reshape(LAB_array,sizey,sizex,3);
XYZ = reshape(XYZ_array,sizey,sizex,3);


% Save the loaded data into the created directories with proper names
save([path '/' Sample_Name '/truth/220 transmittance/transmittance'],'trans_array_m');
save([path '/' Sample_Name '/truth/420 CIEXYZ/XYZ'],'XYZ','sizex','sizey');
save([path '/' Sample_Name '/truth/520 CIELAB/LAB'],'LAB','sizex','sizey');
save([path '/' Sample_Name '/truth/900 sRGB/sRGB'],'rgb');
imwrite(truth,[path '/' Sample_Name '/truth/900 sRGB/truth.tif']);

end