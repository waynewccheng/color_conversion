
classdef MultispectralClass
    
    %% Multispectral Imaging System processing
    % 11-24-2018 Thanksgiving
    
    properties
        Property1
    end
    
    methods (Static)
        
        function obj = MultispectralClass
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
        end
        
        %
        % 11-24-2018 into class
        % cannot find white folder
        % 7-30-2015
        % convert frames (DDL) to reflectance by using reference white background
        %
        function [transmittance_array, sizey, sizex] = frame2reflectance_white (foldername, foldername_white)
            
            disp('Combining frames into reflectance...')
            
            fnin = sprintf('%s/frame',foldername_white);
            load(fnin,'vimarray');
            vimarray0 = vimarray;
            
            [sizewl sizey sizex] = size(vimarray0);
            
            fnin = sprintf('%s/frame',foldername);
            load(fnin,'vimarray');
            
            % calculate the reflectance
            ddl_array = reshape(vimarray,sizewl,sizey*sizex);
            ddl_white_array = reshape(vimarray0,sizewl,sizey*sizex);
            
            transmittance_array = ddl_array ./ ddl_white_array;
            
            transmittance_array = min(transmittance_array,1);
            
            % ----------------------------------
            return
        end
        
        %
        % 11-24-2018
        %
        function spd_array = transmittance2spd (transmittance_array, ls)
            
            vector_length = size(transmittance_array,2);
            ls_array = repmat(ls,1,vector_length);
            spd_array = transmittance_array .* ls_array;
            
        end
        
        
        %
        % Convert the old "reflectance.mat" to "transmittance.mat"
        %
        function convert_transmittance
            
            load('reflectance')
            sizex = 3376
            sizey = 2704
            transmittance_array = reshape(reflectance_array,41,sizey,sizex);
            
            save('transmittance','transmittance_array','sizex','sizey','-v7.3')
            
        end
        
        %
        % Convert the old "reflectance.mat" to "transmittance.mat"
        %
        function convert_transmittance_small
            
            load('reflectance')
            sizex = 844
            sizey = 676
            transmittance_array = reshape(reflectance_array,41,sizey,sizex);
            
            save('transmittance','transmittance_array','sizex','sizey','-v7.3')
            
        end
        
        %
        % Generate monochrome from LAB
        %
        function LAB_mono = LAB2monochrome (LAB)
            LAB_mono = LAB;               % copy LAB
            LAB_mono(:,:,2) = 0;          % set a* = 0
            LAB_mono(:,:,3) = 0;          % set b* = 0
        end
        
        %
        %
        %
        function workflow_dE (tissue_name, pathname_root)
            pathname_tissue = [pathname_root '\' tissue_name];
            pathname_truth = [pathname_tissue '\truth'];
            pathname_mono = [pathname_tissue '\mono'];
            
            pathname_dE = [pathname_tissue '\dE'];
            if exist(pathname_dE,'dir') ~= 7
                mkdir(pathname_dE);
            end
            
            %% LAB_truth in
            folderin = [pathname_truth '\520 CIELAB'];
            load([folderin '\LAB.mat'],'LAB','sizex','sizey')
            LAB_truth = LAB;
            LAB_truth_1d = reshape(LAB_truth,sizey*sizex,3);
            
            %% LAB_mono in
            folderin = [pathname_mono '\500 CIELAB'];
            load([folderin '\LAB.mat'],'LAB_mono','sizex','sizey')
            LAB_mono_1d = reshape(LAB_mono,sizey*sizex,3);
            
            %% dE
            [dE00_1d dE94_1d dEab_1d] = ColorConversionClass.LAB2dE(LAB_truth_1d',LAB_mono_1d');
            dE00 = reshape(dE00_1d,sizey,sizex);
            dE94 = reshape(dE94_1d,sizey,sizex);
            dEab = reshape(dEab_1d,sizey,sizex);
            
            if 0
                %% verify with Sharma's code
                n = sizex * sizey;
                dE00_single = zeros(1,n);
                dE94_single = zeros(1,n);
                dEab_single = zeros(1,n);
                for i = 1:n 
                    [dE00_single(i) dE94_single(i) dEab_single(i)] = ColorConversionClass.LAB2dE(LAB_truth_1d(i,:)',LAB_mono_1d(i,:)');
                    dE00_sharma_single(i) = ColorConversionClass.deltaE2000_Sharma(LAB_truth_1d(i,:),LAB_mono_1d(i,:));
                end
                
                boxplot([dE00_sharma_single' dE00_single' dE94_single' dEab_single' dE00_1d' dE94_1d' dEab_1d'],...
                    {'dE00 Sharma','dE00 single','dE94 single','dEab single','dE00 vector','dE94 vector','dEab vector'})
                
            end
            
            
            %% dE out
            save([pathname_dE '\dE.mat'],'dEab','dE94','dE00','sizex','sizey','-v7.3')
            
        end
        
        
        %
        %
        %
        function workflow_monochrome (tissue_name, pathname_root)
            
            % tissue_name = 'colon';
            % pathname_root = ['C:\Users\wcc\Documents\GitHub\wsi_color_truthing_data'];
            
            pathname_tissue = [pathname_root '\' tissue_name];
            pathname_truth = [pathname_tissue '\truth'];
            
            pathname_mono = [pathname_tissue '\mono'];
            if exist(pathname_mono,'dir') ~= 7
                mkdir(pathname_mono);
            end
            
            %% LAB in
            folderin = [pathname_truth '\520 CIELAB'];
            load([folderin '\LAB.mat'],'LAB','sizex','sizey')
            
            LAB_mono = MultispectralClass.LAB2monochrome(LAB);
            
            %% LAB out
            folderout = [pathname_mono '\500 CIELAB'];
            if exist(folderout,'dir') ~= 7
                mkdir(folderout);
            end
            save([folderout '\LAB.mat'],'LAB_mono','sizex','sizey','-v7.3')
            
            %% sRGB out
            folderout = [pathname_mono '\500 CIELAB'];
            load([folderout '\LAB.mat'],'LAB_mono','sizex','sizey')
            
            %
            % use Matlab lab2rgb() instead of my XYZ2sRGB
            %
            sRGB = lab2rgb(LAB_mono,'ColorSpace','srgb','WhitePoint','d65');
            
            folderout = [pathname_mono '\900 sRGB'];
            if exist(folderout,'dir') ~= 7
                mkdir(folderout);
            end
            
            save([folderout '\sRGB.mat'],'sRGB','sizex','sizey','-v7.3')
            imwrite(sRGB,[folderout '\monochrome.tif'])
            
        end
        
        %
        %
        %
        function workflow_truth (tissue_name, pathname_root)
            
            % tissue_name = 'colon';
            % pathname_root = ['C:\Users\wcc\Documents\GitHub\wsi_color_truthing_data'];
            
            
            pathname_tissue = [pathname_root '\' tissue_name];
            pathname_truth = [pathname_tissue '\truth'];
            
            if 1
                %% light source in
                ls = ColorConversionClass.spd_d65;
                folderout = [pathname_truth '\200 d65'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\d65.SPD.mat'],'ls')
            end
            
            if 1
                %% light source out
                XYZn = ColorConversionClass.spd2XYZ(ColorConversionClass.spd_d65());
                folderout = [pathname_truth '\400 CIEXYZ d65'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\d65.XYZ.mat'],'XYZn')
            end
            
            if 1
                %% transmittance in
                ls = ColorConversionClass.spd_d65();
                folderin = [pathname_truth '\220 transmittance'];
                load([folderin '\transmittance.mat'],'transmittance_array','sizex','sizey')
                transmittance_array_1d = reshape(transmittance_array,41,sizey*sizex);
            end
            
            if 1
                %% SPD out
                spd_array_1d = MultispectralClass.transmittance2spd(transmittance_array_1d, ls);
                spd_array = reshape(spd_array_1d,41,sizey,sizex);
                
                folderout = [pathname_truth '\320 SPD'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\SPD.mat'],'spd_array','sizex','sizey','-v7.3')
            end
            
            if 1
                %% XYZ out
                XYZ_1d = ColorConversionClass.spd2XYZ(spd_array_1d);
                XYZ = reshape(XYZ_1d,sizey,sizex,3);
                
                folderout = [pathname_truth '\420 CIEXYZ'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\XYZ.mat'],'XYZ','sizex','sizey','-v7.3')
            end
            
            if 1
                %% LAB out
                LAB_1d = ColorConversionClass.XYZ2lab(XYZ_1d,XYZn);
                LAB = reshape(LAB_1d,sizey,sizex,3);
                
                folderout = [pathname_truth '\520 CIELAB'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\LAB.mat'],'LAB','sizex','sizey','-v7.3')
            end
            
            if 1
                %% sRGB out
                folderin = [pathname_truth '\420 CIEXYZ'];
                load([folderin '\XYZ.mat'],'XYZ','sizex','sizey')
                XYZ_1d = reshape(XYZ,sizey*sizex,3);
                sRGB_1d = ColorConversionClass.XYZ2sRGB(XYZ_1d);
                sRGB = reshape(sRGB_1d,sizey,sizex,3);
                
                folderout = [pathname_truth '\900 sRGB'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\sRGB.mat'],'sRGB','sizex','sizey','-v7.3')
                
                imwrite(sRGB,[folderout '\truth.tif'])
            end
            
            
        end
        
    end
end

