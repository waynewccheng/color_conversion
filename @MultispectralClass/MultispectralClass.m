
classdef MultispectralClass
    
    %% Multispectral Imaging System processing
    % 11-24-2018 Thanksgiving
    
    properties
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
        % convert frames (DDL) to transmittance by using reference white background
        %
        function [transmittance_array, sizey, sizex] = frame2transmittance_white_black (vimarray, vimarray_white, vimarray_black)
            
            disp('Combining frames into transmittance...')
            
            [sizewl sizey sizex] = size(vimarray_white);
            
            % calculate the transmittance
            ddl_array = reshape(vimarray,sizewl,sizey*sizex);
            ddl_white_array = reshape(vimarray_white,sizewl,sizey*sizex);
            ddl_black_array = reshape(vimarray_black,sizewl,sizey*sizex);
            
if 0            
            % set the transmittance as 0 if the reference white is too dark
            threshold = 300;
            
             mask = (ddl_white_array - ddl_black_array) <= threshold |...
             ddl_white_array < threshold |...
             ddl_array < ddl_black_array;
                 
             transmittance_array(mask) = 0;
             transmittance_array(~mask) = ddl_array(~mask) ./ ddl_white_array(~mask);
    end
    
if 1
            transmittance_array = (ddl_array - ddl_black_array) ./ (ddl_white_array - ddl_black_array);

            threshold = 150;
            
            mask_greater_than_white = (ddl_array >= ddl_white_array) & ((ddl_white_array - ddl_black_array) > threshold);
            transmittance_array(mask_greater_than_white) = 1;

            mask_less_than_black = (ddl_array <= ddl_black_array) & ((ddl_white_array - ddl_black_array) > threshold);
            transmittance_array(mask_less_than_black) = 0;

            mask_dynamic_range_too_small = (ddl_white_array - ddl_black_array) < threshold;
            transmittance_array(mask_dynamic_range_too_small) = 0;
            
end
if 0
            % ----------------------------------
            % clip at 1 => not a good solution !!!
            transmittance_array = (ddl_array - ddl_black_array) ./ (ddl_white_array - ddl_black_array);
            transmittance_array = min(transmittance_array,1);
            transmittance_array = max(transmittance_array,0);
end


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
                save([pathname_truth '\info.mat'],'xy','z','xy_white','z_white','xy_black','z_black')
                
            else
                
                % shortcut: use the previously saved locations
                load([pathname_truth '\info.mat'],'xy','z','xy_white','z_white','xy_black','z_black')
                
            end
            
            %% Step 1: take frames
            
            if 1
                
                n_times = 10;
                
                disp('Moving to the ROI')
                ludl.setXY(xy)
                ludl.setZ(z)
                
                disp('Getting frames for ROI')
                vimarray = MultispectralClass.camera2frame(pathname_truth,n_times,ol490,big_or_small);
                
                disp('Moving to the reference white')
                ludl.setXY(xy_white)
                ludl.setZ(z_white)
                
                disp('Getting frames for white')
                vimarray_white = MultispectralClass.camera2frame(pathname_truth,n_times,ol490,big_or_small);

                disp('Moving to the reference black')
                ludl.setXY(xy_black)
                ludl.setZ(z_black)
                
                disp('Getting frames for black')
                vimarray_black = MultispectralClass.camera2frame(pathname_truth,n_times,ol490,big_or_small);
                
                
                ol490.setWhite      % recover the light to white
                
                ludl.setXY(xy)      % return to the ROI
                ludl.setZ(z)
                
                ludl.close;
                
                
                % save data after closing devices
                disp('Saving captured frames in vimarray...')
                
                tic
                folderout = [pathname_truth '\120 frame'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                
                var_size = whos('vimarray');
                if var_size.bytes > 2*1000*1000*1000
                    save([folderout '\frame.mat'],'vimarray','-v7.3','-nocompression')
                else
                    save([folderout '\frame.mat'],'vimarray')
                end
                toc
                
                tic
                folderout = [pathname_truth '\110 frame white'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                
                var_size = whos('vimarray_white');
                if var_size.bytes > 2*1000*1000*1000
                    save([folderout '\white.frame.mat'],'vimarray_white','-v7.3','-nocompression')
                else
                    save([folderout '\white.frame.mat'],'vimarray_white')
                end
                toc

                tic
                folderout = [pathname_truth '\100 frame black'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                
                var_size = whos('vimarray_black');
                if var_size.bytes > 2*1000*1000*1000
                    save([folderout '\black.frame.mat'],'vimarray_black','-v7.3','-nocompression')
                else
                    save([folderout '\black.frame.mat'],'vimarray_black')
                end
                toc                
                
                return
                
            end
            
            %
            % 4-10-2018
            %
            function findroi
                
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
                
                a = input('Press Enter to save location of black:')
                xy_black = ludl.getXY
                z_black = ludl.getZ
                
                delete(vid)
                
            end
            
        end
        
        %% Show images of 41 bands in subplots
        %
        % Q: Is the light field uniform?
        %
        % filename: .mat, use fullname including path and extension
        % varname: vimarray or vimarray_white
        %
        % revised 7-18-2018: added under- and over-flow pixels in percent
        % revised 11-27-2018: into class
        %
        function frame2plot (foldername, filename, varname)
            
            load([foldername '\' filename],varname)
            
            sizey = size(eval(varname),2);
            sizex = size(eval(varname),3);
            
            % the "reflectance" matrix is 41x570544
            whos
            
            % get the pixel count
            n_total = sizey * sizex;
            
            % convert the matrix from 1D to 2D
            data = reshape(eval(varname),41,sizey,sizex);
            
            % create a big figure; otherwise the titles are unreadable
            h1 = figure('position',[10 100 1900 1000]);
            
            % sweep the wavelength
            wl = 380;
            for i = 1:41
                
                % retrieve the frame
                im1 = squeeze(data(i,:,:));
                
                % get statistics
                n_over =  nnz(im1 >= 1);
                n_under = nnz(im1 <= 0);
                
                % re-create a color image
                %         im = uint8(zeros(sizey,sizex,3));
                %         im(:,:,1) = im1*255;
                %         im(:,:,2) = im1*255;
                %         im(:,:,3) = im1*255;
                
                % find a space to show the image
                subplot(6,7,i)
                
                imagesc(im1)
                axis off
                
                % show the statistics
                title(sprintf('%d: [%.2f%% %.2f%%]',wl,n_under/n_total*100,n_over/n_total*100),'FontSize',8)
                
                % post iteration
                wl = wl + 10;
            end
            
            % save the image
            saveas(h1,[foldername '\frames41.png'])
            savefig(h1,[foldername '\frames41'])  % the png is awful
        end
        
        
        %% Show multispectral images of 41 bands
        % revised 7-18-2018: added under- and over-flow pixels in percent
        %
        function frame2boxplot (foldername, filename, varname)
            
            filepath = [foldername '\' filename];
            load(filepath,varname)
            
            sizey = size(eval(varname),2);
            sizex = size(eval(varname),3);
            n_total = sizey * sizex;
            
            vimarray_1d = reshape(eval(varname),41,n_total);
            
            % the "reflectance" matrix is 41x570544
            whos
            
            % create a big figure; otherwise the titles are unreadable
            h1 = figure('position',[10 100 1900 1000]);
            
            % boxplot
            boxplot(vimarray_1d')
            
            % save the image
            saveas(h1,[foldername '\boxplot.png'])
        end
        
        
        %
        %
        %
        function transmittance2boxplot (pathname_truth)
            
            % read the "reflectance.mat" file from the folder
            load([pathname_truth '\220 transmittance\transmittance'],'transmittance_array','sizey','sizex')
            
            % the "reflectance" matrix is 41x570544
            whos
            
            % create a big figure; otherwise the titles are unreadable
            h1 = figure('position',[10 100 1900 1000]);
            
            % boxplot
            boxplot(transmittance_array')
            
            % save the image
            % save the image
            folderout = [pathname_truth '\220 transmittance'];
            if exist(folderout,'dir') ~= 7
                mkdir(folderout);
            end
            saveas(h1,[folderout '\boxplot.png'])
            
        end
        
        %% Visually check the transmittance and frames
        %
        %
        function workflow_acquire_sanity_check (pathname_truth)
            
            if 0
                MultispectralClass.frame2boxplot([pathname_truth '\120 frame'],'frame.mat','vimarray')
                MultispectralClass.frame2boxplot([pathname_truth '\110 frame white'],'white.frame.mat','vimarray_white')
                MultispectralClass.frame2boxplot([pathname_truth '\100 frame black'],'black.frame.mat','vimarray_black')
                MultispectralClass.frame2plot([pathname_truth '\120 frame'],'frame.mat','vimarray')
                MultispectralClass.frame2plot([pathname_truth '\110 frame white'],'white.frame.mat','vimarray_white')
                MultispectralClass.frame2plot([pathname_truth '\100 frame black'],'black.frame.mat','vimarray_black')
            end
            
            MultispectralClass.transmittance2boxplot(pathname_truth)
            
            disp('Use "close all" to close all figures')
        end
        
        %% Acquire frames and calculate spectral transmittance
        %
        %
        function workflow_acquire (pathname_truth, ol490)
            
            % acquire the frames
            if 0
            MultispectralClass.acquire_frames(pathname_truth,ol490,0)
            end
            
            % calculate the spectral transmittance
            if 1
                %% frame in
                folderin = [pathname_truth '\100 frame black'];
                load([folderin '\black.frame.mat'],'vimarray_black')
                
                folderin = [pathname_truth '\110 frame white'];
                load([folderin '\white.frame.mat'],'vimarray_white')
                
                folderin = [pathname_truth '\120 frame'];
                load([folderin '\frame.mat'],'vimarray')

                [transmittance_array, sizey, sizex] = MultispectralClass.frame2transmittance_white_black(vimarray, vimarray_white, vimarray_black);
                
                folderout = [pathname_truth '\220 transmittance'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '\transmittance.mat'],'transmittance_array','sizex','sizey','-V7.3')
            end
            
        end
        
        %% Calculate SPD, XYZ, LAB, and sRGB
        %
        %
        function workflow_truth (pathname_truth)
            
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
        
        
        %% Calculate monochrome
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
        
        %% Calculate dE
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
        
        
    end
end

