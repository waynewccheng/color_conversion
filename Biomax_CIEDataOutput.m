% Biomax Data Output
% Given transmittance data output from Paul's Code (for Biomax samples),
% this code will output SPD, XYZ, LAB, and sRGB into the proper folders

% Make sure to change sample name (Line 8) each time

path = ('E:/DigitalPathology/HIMS_Data/ProcessedData/Truth/');
sample = 'BiomaxOrgan10_Bladder_M13';
path_truth = [path sample];

%% light source out
        XYZn = ColorConversionClass.spd2XYZ(ColorConversionClass.spd_d65());
        
%% transmittance in
        ls = ColorConversionClass.spd_d65();
        folderin = [path_truth '/Transmittance'];
        load([folderin '/trans_mean_camera.mat'],'trans_array_m','sizex','sizey')
%         transmittance_array_1d = reshape(trans_array_m,41,sizey*sizex);
        
%% SPD out
        spd_array_1d = MultispectralClass.transmittance2spd(trans_array_m, ls);
        spd_array = reshape(spd_array_1d,41,sizey,sizex);   
        folderout = [path_truth '/CIE_Coord'];
        save([folderout '/SPD.mat'],'spd_array','sizex','sizey','-v7.3')
        
%% XYZ out
        XYZ_array = ColorConversionClass.spd2XYZ(spd_array_1d);
        XYZ = reshape(XYZ_array,sizey,sizex,3);              
        folderout = [path_truth '/CIE_Coord'];
        save([folderout '/XYZ_array.mat'],'XYZ_array')

%% LAB out
        LAB_array = ColorConversionClass.XYZ2lab(XYZ_array,XYZn);
        LAB = reshape(LAB_array,sizey,sizex,3);       
        folderout = [path_truth '/CIE_Coord'];
        save([folderout '/LAB_array.mat'],'LAB_array')

%% sRGB out
        rgb = ColorConversionClass.XYZ2sRGB(XYZ_array);
        sRGB = reshape(rgb,sizey,sizex,3);      
        folderout = [path_truth '/EndResults'];
        save([folderout '/rgb.mat'],'rgb')