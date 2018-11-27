
classdef MultispectralClass
    
    %% Multispectral Imaging System processing
    % 11-24-2018 Thanksgiving
    
    properties
        Property1
    end
    
    methods (Static)
        
        function obj = MultispectralClass
        end
        
        %
        % 4-12-2018: save cam.vid and cam.src in info.mat
        % 4-10-2018: merge big and small
        % 10-9-2015
        % example: ol490 = OL490Class; camera2frame('dataout/1009-test1',1,ol490)
        %
        % 8-3-2015: hardware trigger; needs software reset!
        % 7-23-2015: replace grasshopper function with camera_*
        % usage: camera2frame('0723-7')
        % 7-21-2015: revisit
        % capture images with camera
        % capture 41 images from 380 to 780
        % output: vimarray(41,480,640)
        %
        function vimarray = camera2frame (pathout, numberofshots, ol490, big_or_small)
            
            disp('Capturing frames with wavelength from 380 to 780 nm...')
            
            
            % aqusition
            if big_or_small == 1
                cam = CameraClass9MPBig
            else
                cam = CameraClass9MPSmall
            end
            
            cam_vid = get(cam.vid)
            cam_src = get(cam.src)
            
            %            mkdir(pathout)
            %            fnout_info = sprintf('%s/info',pathout);
            %            save(fnout_info,'cam_vid','cam_src','-append')
            
            % data
            vimarray = zeros(41,cam.sizey,cam.sizex);
            
            
            % prepare light
            bandwidth = 10;
            intensity = 100;
            
            % add some delay here because 380 nm has problems
            ol490.setPeak(380,bandwidth,intensity);
            pause(1)
            
            k = 1;
            for wl=380:10:780
                % prepare light
                
                ol490.setPeak(wl,bandwidth,intensity);
                
                % focus
                %        f_opt = myfocus(cam)
                
                % need pause here???
                % pause(0.25)
                
                % acqusition
                vim = cam.snap(numberofshots);
                vimarray(k,:,:) = vim;
                
                k = k + 1;
            end
            
            % exit
            cam.close;
            
            beep
            
            % save data after closing devices
            %disp('Saving captured frames in vimarray...')
            %save(fnout,'vimarray','-V7.3')
            
            ol490.setPeak(550,10,100)
            
            % ring
            beep on
            beep
            
            return
            
        end
        
        %
        % 11-24-2018 into class
        % cannot find white folder
        % 7-30-2015
        % convert frames (DDL) to reflectance by using reference white background
        %
        function [transmittance_array, sizey, sizex] = frame2reflectance_white (vimarray, vimarray0)
            
            disp('Combining frames into reflectance...')
            
            [sizewl sizey sizex] = size(vimarray0);
            
            % calculate the reflectance
            ddl_array = reshape(vimarray,sizewl,sizey*sizex);
            ddl_white_array = reshape(vimarray0,sizewl,sizey*sizex);
            
            transmittance_array = ddl_array ./ ddl_white_array;
            
            % clip at 1
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
        function workflow_monochrome (pathname_truth, pathname_mono)
            
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
        function acquire_frames (pathname_truth, ol490, big_or_small)
            
            tic
            
            % load('config.mat','LUDL_PORT');
            LUDL_PORT = 'COM14'
            
            ludl = LudlClass(LUDL_PORT);
            
            %% Step 0: choose ROI and focus
            
            % 0: use new location; 1: use old location
            if 0
                
                % determine the ROI locations for target and reference white
                findroi()
                
                % save the ROI locations
                if exist(pathname_truth,'dir') ~= 7
                    mkdir(pathname_truth)
                end
                save([pathname_truth '\info.mat'],'xy','xy_white','z','z_white')
                
            else
                
                % shortcut: use the previously saved locations
                load('info.mat','xy','xy_white','z','z_white')
                
            end
            
            %% Step 1: take frames
            
            if 1
                
                disp('Moving to the ROI')
                ludl.setXY(xy)
                ludl.setZ(z)
                
                disp('Getting frames for ROI')
                vimarray = MultispectralClass.camera2frame(pathname_truth,1,ol490,big_or_small);
                
                disp('Moving to the reference white')
                ludl.setXY(xy_white)
                ludl.setZ(z_white)
                
                disp('Getting frames for white')
                vimarray0 = MultispectralClass.camera2frame(pathname_truth,1,ol490,big_or_small);
                
                
                ol490.setWhite      % recover the light to white
                
                ludl.setXY(xy)      % return to the ROI
                ludl.setZ(z)
                
                ludl.close;
                
                
                % save data after closing devices
                disp('Saving captured frames in vimarray...')
                
                folderout = [pathname_truth '\120 frame'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\frame.mat'],'vimarray','-V7.3')
                
                folderout = [pathname_truth '\110 frame white'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\white.frame.mat'],'vimarray0','-V7.3')
                
            end
            
            %
            % 4-10-2018
            
            % turn on the light
            function findroi
                
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
            end
            
        end
        
        %
        %
        %
        function workflow_truth (pathname_truth)
            
            if 1
                %% frame in
                folderin = [pathname_truth '\110 frame white']
                load([folderin '\white.frame.mat'],'vimarray0')
                
                folderin = [pathname_truth '\120 frame'];
                load([folderin '\frame.mat'],'vimarray')
                
                [transmittance_array, sizey, sizex] = MultispectralClass.frame2reflectance_white(vimarray, vimarray0);
                
                folderout = [pathname_truth '\220 transmittance'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\transmittance.mat'],'transmittance_array','sizex','sizey','-V7.3')
            end
            
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

