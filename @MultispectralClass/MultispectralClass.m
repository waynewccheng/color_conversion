
classdef MultispectralClass
    
    %% Multispectral Imaging System processing
    % 11-24-2018 Thanksgiving
    %
    % 09-20-2020 ACL -- Modified version of MultispectralClass to register images
    % from different scanners and compute dE between them
    %
    % This version of the code assumes that, for a given sample, you have
    % transmittance, LAB, XYZ, and RGB data output by Paul's HIMS code and
    % have a WSI image of the same sample from different scanners.
    %
    % General Workflow for a given Scanner_Name:
    % 1) Modify definitions in conversion function (Lines 1029-1041) as needed
    % 2) Run MultispectralClass.conversion in command window to create folders & convert files 
    % 3) Manually save an image of the ROI from the scanner as "scan.tif" in the "Scanner_Name/200 roi" folder location
    % 4) Use the Registration Estimator App to register your truth image with your ROI image from the scanner
    % 5) Save the output 'movingReg' file as 'matlab_registration' in the "Scanner_Name/400 sRGB" folder
    % 6) Modify definitions in register_scanner_tissue function (Lines 1013-1024) as needed
    % 7) Run MultispectralClass.register_scanner_tissue in command window
    
    properties
    end
    
    methods (Static)
        %register_controlpoint_save (original_fn, unregistered_fn, folder_output)
        [A,B] = register_controlpoint_show (original_fn, unregistered_fn, folder_output, lab_truth_fn)
        check_registration_with_bars (folder)
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
        function [vimarray caminfo] = camera2frame (numberofshots, ol490, big_or_small)
            
            disp('Capturing frames with wavelength from 380 to 780 nm...')
            
            
            % aqusition
            if big_or_small == 1
                cam = CameraClass9MPBig
            else
                cam = CameraClass9MPSmall
            end
            
            caminfo.vid = get(cam.vid);
            caminfo.src = get(cam.src);
            caminfo.sizex = cam.sizex;
            caminfo.sizey = cam.sizey;
            caminfo.time1 = datetime;
            
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
                
                fprintf('%d ',wl)
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
            
            caminfo.time2 = datetime;
            disp(' done!')
            
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
            
            % prepare the ROI, white, and black frames
            ddl_array = reshape(vimarray,sizewl,sizey*sizex);
            ddl_white_array = reshape(vimarray_white,sizewl,sizey*sizex);
            ddl_black_array = reshape(vimarray_black,sizewl,sizey*sizex);
            
            % analyze the frame
            analyze_frame('order_stat_before');
            
            % condition the frames
            
            %  analyze the frame
            analyze_frame('order_stat_after');
            
            % calculate the transmittance
            transmittance_array = (ddl_array - ddl_black_array) ./ (ddl_white_array - ddl_black_array);
            
            
            %% Condition the transmittance
            if 1
                % TBD
                threshold = 150;
                
                mask_greater_than_white = (ddl_array >= ddl_white_array) & ((ddl_white_array - ddl_black_array) > threshold);
                transmittance_array(mask_greater_than_white) = 1;
                
                mask_less_than_black = (ddl_array <= ddl_black_array) & ((ddl_white_array - ddl_black_array) > threshold);
                transmittance_array(mask_less_than_black) = 0;
                
                mask_dynamic_range_too_small = (ddl_white_array - ddl_black_array) < threshold;
                transmittance_array(mask_dynamic_range_too_small) = 0;
                
            else
                % ----------------------------------
                % clip at 1 => not a good solution !!!
                transmittance_array = (ddl_array - ddl_black_array) ./ (ddl_white_array - ddl_black_array);
                transmittance_array = min(transmittance_array,1);
                transmittance_array = max(transmittance_array,0);
            end
            
            return
            
            %
            %
            %
            function analyze_frame (fn)
                
                if 1
                    mask_bxw = ddl_black_array <= ddl_array & ddl_array <= ddl_white_array;
                    mask_bwx = ddl_black_array <= ddl_white_array & ddl_white_array <= ddl_array;
                    mask_xbw = ddl_array <= ddl_black_array & ddl_black_array <= ddl_white_array;
                    mask_xwb = ddl_array <= ddl_white_array & ddl_white_array <= ddl_black_array;
                    mask_wbx = ddl_white_array <= ddl_black_array & ddl_black_array <= ddl_array;
                    mask_wxb = ddl_white_array <= ddl_array & ddl_array <= ddl_black_array;
                else
                    mask_bxw = ddl_black_array < ddl_array & ddl_array < ddl_white_array;
                    mask_bwx = ddl_black_array < ddl_white_array & ddl_white_array < ddl_array;
                    mask_xbw = ddl_array < ddl_black_array & ddl_black_array < ddl_white_array;
                    mask_xwb = ddl_array < ddl_white_array & ddl_white_array < ddl_black_array;
                    mask_wbx = ddl_white_array < ddl_black_array & ddl_black_array < ddl_array;
                    mask_wxb = ddl_white_array < ddl_array & ddl_array < ddl_black_array;
                end
                
                order_stat = zeros(6,41);
                for i = 1:41
                    order_stat(1,i) = nnz(mask_bxw(i,:));
                    order_stat(2,i) = nnz(mask_bwx(i,:));
                    order_stat(3,i) = nnz(mask_xbw(i,:));
                    order_stat(4,i) = nnz(mask_xwb(i,:));
                    order_stat(5,i) = nnz(mask_wbx(i,:));
                    order_stat(6,i) = nnz(mask_wxb(i,:));
                end
                
                order_stat = order_stat / size(ddl_array,2);
                
                save(fn,'order_stat')
                
                % visualization
                clf
                subplot(2,2,4)
                bar(order_stat','stacked')
                legend('bxw','bwx','xbw','xwb','wbx','wxb','Location','best')
                
                subplot(2,2,1)
                hold on
                plot(order_stat(1,:),'r-')
                plot(order_stat(2,:),'r:')
                plot(order_stat(3,:),'g-')
                plot(order_stat(4,:),'g:')
                plot(order_stat(5,:),'b-')
                plot(order_stat(6,:),'b:')
                axis([1 41 0 1])
                legend('bxw','bwx','xbw','xwb','wbx','wxb','Location','best')
                
                subplot(2,2,2)
                hold on
                plot(order_stat(1,:),'r-')
                plot(order_stat(2,:),'r:')
                plot(order_stat(3,:),'g-')
                plot(order_stat(4,:),'g:')
                plot(order_stat(5,:),'b-')
                plot(order_stat(6,:),'b:')
                axis([1 3 0 0.01])
                legend('bxw','bwx','xbw','xwb','wbx','wxb','Location','best')
                
                subplot(2,2,3)
                hold on
                plot(order_stat(1,:),'r-')
                plot(order_stat(2,:),'r:')
                plot(order_stat(3,:),'g-')
                plot(order_stat(4,:),'g:')
                plot(order_stat(5,:),'b-')
                plot(order_stat(6,:),'b:')
                axis([41-2 41 0 0.01])
                legend('bxw','bwx','xbw','xwb','wbx','wxb','Location','best')
                
            end
            
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
        function acquire_frames (location_filename, pathname_truth, ol490, big_or_small)
            
            % load('config.mat','LUDL_PORT');
            LUDL_PORT = 'COM14'
            ludl = LudlClass(LUDL_PORT);
            
            %% Step 0: choose ROI and focus
            
            % bypass scanning: use the previously saved locations
            if exist(location_filename,'file')
                load(location_filename,'xyz*')
            else
                disp('Location file not found!!!!!')
                exit
            end
            
            %% Step 1: take frames
            
            if 1
                
                n_times = 10;
                
                disp('Moving to the ROI')
                tic
                ludl.setXY(xyz(1:2))
                ludl.setZ(xyz(3))
                toc
                
                disp('Getting frames for ROI')
                tic
                [vimarray caminfo] = MultispectralClass.camera2frame(n_times,ol490,big_or_small);
                toc
                
                disp('Moving to the reference white')
                tic
                ludl.setXY(xyz_white(1:2))
                ludl.setZ(xyz_white(3))
                toc
                
                disp('Getting frames for white')
                tic
                [vimarray_white caminfo_white] = MultispectralClass.camera2frame(n_times,ol490,big_or_small);
                toc
                
                disp('Moving to the reference black')
                tic
                ludl.setXY(xyz_black(1:2))
                ludl.setZ(xyz_black(3))
                toc
                
                disp('Getting frames for black')
                tic
                [vimarray_black caminfo_black] = MultispectralClass.camera2frame(n_times,ol490,big_or_small);
                toc
                
                ol490.setWhite      % recover the light to white
                
                ludl.setXY(xyz(1:2))      % return to the ROI
                ludl.setZ(xyz(3))
                
                ludl.close;
                
                
                % save data after closing devices
                disp('Saving truth frames')
                tic
                folderout = [pathname_truth '/120 frame'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                
                var_size = whos('vimarray');
                if var_size.bytes > 2*1000*1000*1000
                    %                    save([folderout '/frame.mat'],'vimarray','-v7.3','-nocompression')
                    save([folderout '/frame.mat'],'vimarray','-v7.3')
                else
                    save([folderout '/frame.mat'],'vimarray')
                end
                save([folderout '/caminfo.mat'],'caminfo','xyz')
                toc
                
                disp('Saving white frames')
                tic
                folderout = [pathname_truth '/110 frame white'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                
                var_size = whos('vimarray_white');
                if var_size.bytes > 2*1000*1000*1000
                    %                    save([folderout '/white.frame.mat'],'vimarray_white','-v7.3','-nocompression')
                    save([folderout '/white.frame.mat'],'vimarray_white','-v7.3')
                else
                    save([folderout '/white.frame.mat'],'vimarray_white')
                end
                save([folderout '/caminfo.mat'],'caminfo_white','xyz_white')
                toc
                
                disp('Saving black frames')
                tic
                folderout = [pathname_truth '/100 frame black'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                
                var_size = whos('vimarray_black');
                if var_size.bytes > 2*1000*1000*1000
                    %                    save([folderout '/black.frame.mat'],'vimarray_black','-v7.3','-nocompression')
                    save([folderout '/black.frame.mat'],'vimarray_black','-v7.3')
                else
                    save([folderout '/black.frame.mat'],'vimarray_black')
                end
                save([folderout '/caminfo.mat'],'caminfo_black','xyz_black')
                toc
                
                return
                
            end
            
            
        end

        % define ROI locations
        % 11-30-2018
        % 4-10-2018
        function find_ROI_white_black (ol490, ludl_port_name, filename)
            
            %ludl_port_name = 'COM14'
            
            % turn on the light
            ol490.setWhite
            
            vid = videoinput('pointgrey', 1, 'F7_Mono8_844x676_Mode5');
            src = getselectedsource(vid);
            
            vid.FramesPerTrigger = 1;
            
            preview(vid);
            
            ludl = LudlClass(ludl_port_name)
            
            
            a = input('Press Enter to save location of ROI:')
            
            xy = ludl.getXY
            z = ludl.getZ
            xyz = [xy z]
            
            
            a = input('Press Enter to save location of white:')
            
            xy_white = ludl.getXY
            z_white = ludl.getZ
            xyz_white = [xy_white z_white]
            
            
            a = input('Press Enter to save location of black:')
            
            xy_black = ludl.getXY
            z_black = ludl.getZ
            xyz_black = [xy_black z_black]
            
            ludl.close
            
            delete(vid)
            
            save(filename,'xyz*')
            
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
            
            load([foldername '/' filename],varname)
            
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
            saveas(h1,[foldername '/frames41.png'])
%            savefig(h1,[foldername '/frames41'])  % the png is awful
        end
        
        
        %% Show multispectral images of 41 bands
        % revised 7-18-2018: added under- and over-flow pixels in percent
        %
        function frame2boxplot (foldername, filename, varname)
            
            filepath = [foldername '/' filename];
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
            saveas(h1,[foldername '/boxplot.png'])
        end
        
        
        %
        %
        %
        function transmittance2boxplot (pathname_truth)
            
            % read the "reflectance.mat" file from the folder
            load([pathname_truth '/220 transmittance/transmittance'],'transmittance_array','sizey','sizex')
            
            % the "reflectance" matrix is 41x570544
            whos
            
            % create a big figure; otherwise the titles are unreadable
            h1 = figure('position',[10 100 1900 1000]);
            
            % boxplot
            boxplot(transmittance_array')
            
            % save the image
            % save the image
            folderout = [pathname_truth '/220 transmittance'];
            if exist(folderout,'dir') ~= 7
                mkdir(folderout);
            end
            saveas(h1,[folderout '/boxplot.png'])
            
        end
        
        %% Visually check the transmittance and frames
        %
        %
        function workflow_acquire_sanity_check (pathname_truth)
            
            if 1
                MultispectralClass.frame2boxplot([pathname_truth '/120 frame'],'frame.mat','vimarray')
                MultispectralClass.frame2boxplot([pathname_truth '/110 frame white'],'white.frame.mat','vimarray_white')
                MultispectralClass.frame2boxplot([pathname_truth '/100 frame black'],'black.frame.mat','vimarray_black')
                MultispectralClass.frame2plot([pathname_truth '/120 frame'],'frame.mat','vimarray')
                MultispectralClass.frame2plot([pathname_truth '/110 frame white'],'white.frame.mat','vimarray_white')
                MultispectralClass.frame2plot([pathname_truth '/100 frame black'],'black.frame.mat','vimarray_black')
            end
            
            MultispectralClass.transmittance2boxplot(pathname_truth)
            
            disp('Use "close all" to close all figures')
        end
        
        %% Acquire frames and calculate spectral transmittance
        %
        %
        function workflow_acquire (location_filename, pathname_truth, ol490)
            
            % acquire the frames
            if 1
                % for testing
                % MultispectralClass.acquire_frames(location_filename, pathname_truth,ol490,0) %small
                
                % for testing
                MultispectralClass.acquire_frames(location_filename, pathname_truth,ol490,1)   %big
            end
            
            % calculate the spectral transmittance
            if 1
                %% frame in
                folderin = [pathname_truth '/100 frame black'];
                load([folderin '/black.frame.mat'],'vimarray_black')
                
                folderin = [pathname_truth '/110 frame white'];
                load([folderin '/white.frame.mat'],'vimarray_white')
                
                folderin = [pathname_truth '/120 frame'];
                load([folderin '/frame.mat'],'vimarray')
                
                [transmittance_array, sizey, sizex] = MultispectralClass.frame2transmittance_white_black(vimarray, vimarray_white, vimarray_black);
                
                folderout = [pathname_truth '/220 transmittance'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '/transmittance.mat'],'transmittance_array','sizex','sizey','-V7.3')
            end
            
        end
        
        %% Calculate SPD, XYZ, LAB, and sRGB
        %
        %
        function workflow_truth (pathname_truth)
            
            if exist(pathname_truth,'dir') ~= 7
                mkdir(pathname_truth);
            end
            
            if 1
                %% light source in
                ls = ColorConversionClass.spd_d65;
                folderout = [pathname_truth '/200 d65'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '/d65.SPD.mat'],'ls')
            end
            
            if 1
                %% light source out
                XYZn = ColorConversionClass.spd2XYZ(ColorConversionClass.spd_d65());
                folderout = [pathname_truth '/400 CIEXYZ d65'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '/d65.XYZ.mat'],'XYZn')
            end
            
            if 1
                %% transmittance in
                ls = ColorConversionClass.spd_d65();
                folderin = [pathname_truth '/220 transmittance'];
                load([folderin '/transmittance.mat'],'transmittance_array','sizex','sizey')
                transmittance_array_1d = reshape(transmittance_array,41,sizey*sizex);
            end
            
            if 1
                %% SPD out
                spd_array_1d = MultispectralClass.transmittance2spd(transmittance_array_1d, ls);
                spd_array = reshape(spd_array_1d,41,sizey,sizex);
                
                folderout = [pathname_truth '/320 SPD'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '/SPD.mat'],'spd_array','sizex','sizey','-v7.3')
            end
            
            if 1
                %% XYZ out
                XYZ_1d = ColorConversionClass.spd2XYZ(spd_array_1d);
                XYZ = reshape(XYZ_1d,sizey,sizex,3);
                
                folderout = [pathname_truth '/420 CIEXYZ'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '/XYZ.mat'],'XYZ','sizex','sizey','-v7.3')
            end
            
            if 1
                %% LAB out
                LAB_1d = ColorConversionClass.XYZ2lab(XYZ_1d,XYZn);
                LAB = reshape(LAB_1d,sizey,sizex,3);
                
                folderout = [pathname_truth '/520 CIELAB'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '/LAB.mat'],'LAB','sizex','sizey','-v7.3')
            end
            
            if 1
                %% sRGB out
                folderin = [pathname_truth '/420 CIEXYZ'];
                load([folderin '/XYZ.mat'],'XYZ','sizex','sizey')
                XYZ_1d = reshape(XYZ,sizey*sizex,3);
                sRGB_1d = ColorConversionClass.XYZ2sRGB(XYZ_1d);
                sRGB = reshape(sRGB_1d,sizey,sizex,3);
                
                folderout = [pathname_truth '/900 sRGB'];
                if exist(folderout,'dir') ~= 7
                    mkdir(folderout);
                end
                save([folderout '/sRGB.mat'],'sRGB','sizex','sizey','-v7.3')
                
                imwrite(sRGB,[folderout '/truth.tif'])
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
            folderin = [pathname_truth '/520 CIELAB'];
            load([folderin '/LAB.mat'],'LAB','sizex','sizey')
            
            LAB_mono = MultispectralClass.LAB2monochrome(LAB);
            
            %% LAB out
            folderout = [pathname_mono '/500 CIELAB'];
            if exist(folderout,'dir') ~= 7
                mkdir(folderout);
            end
            save([folderout '/LAB.mat'],'LAB_mono','sizex','sizey','-v7.3')
            
            %% sRGB out
            folderout = [pathname_mono '/500 CIELAB'];
            load([folderout '/LAB.mat'],'LAB_mono','sizex','sizey')
            
            %
            % use Matlab lab2rgb() instead of my XYZ2sRGB
            %
            sRGB = lab2rgb(LAB_mono,'ColorSpace','srgb','WhitePoint','d65');
            
            folderout = [pathname_mono '/900 sRGB'];
            if exist(folderout,'dir') ~= 7
                mkdir(folderout);
            end
            
            save([folderout '/sRGB.mat'],'sRGB','sizex','sizey','-v7.3')
            imwrite(sRGB,[folderout '/monochrome.tif'])
            
        end
        
        %% Calculate dE
        %
        %
        function workflow_dE (pathname_truth, pathname_mono, pathname_dE)
            
            if exist(pathname_dE,'dir') ~= 7
                mkdir(pathname_dE);
            end
            
            %% LAB_truth in
            folderin = [pathname_truth '/520 CIELAB'];
            load([folderin '/LAB.mat'],'LAB','sizex','sizey')
            LAB_truth = LAB;
            LAB_truth_1d = reshape(LAB_truth,sizey*sizex,3);
            
            %% LAB_mono in
            folderin = [pathname_mono '/500 CIELAB'];
            load([folderin '/LAB.mat'],'LAB_mono','sizex','sizey')
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
            
            %% Plot dE heatmap
            %Added 9_30_20 by ACL
            
            X = 1:sizex;
            Y = 1:sizey;
            figure;
            surf(X,Y,dE00,'FaceColor','interp','EdgeColor','none')
            title('\DeltaE00 Monochromatic');
            zlim([0 100]); xlim([0 sizex]); ylim([0 sizey]);
            
        end

        %% Calculate dE for scanner
        %
        %
        function workflow_dE_scan (LAB_truth, LAB_scan, pathname_dE, scanner)
            
            if exist(pathname_dE,'dir') ~= 7
                mkdir(pathname_dE);
            end
            
            %% LAB_truth in
            sizex = size(LAB_truth,2);
            sizey = size(LAB_truth,1);
            LAB_truth_1d = reshape(LAB_truth,sizey*sizex,3);
            
            %% LAB_scan in
            sizex = size(LAB_scan,2);
            sizey = size(LAB_scan,1);
            LAB_scan_1d = reshape(LAB_scan,sizey*sizex,3);
            
            %% dE
            [dE00_1d dE94_1d dEab_1d] = ColorConversionClass.LAB2dE(LAB_truth_1d',LAB_scan_1d');
            dE00 = reshape(dE00_1d,sizey,sizex);
            dE94 = reshape(dE94_1d,sizey,sizex);
            dEab = reshape(dEab_1d,sizey,sizex);
            
            
            %% dE out
            save([pathname_dE '\dE_' scanner '.mat'],'dEab','dE94','dE00','sizex','sizey','-v7.3')
            
            %% plot heatmap of dE
            % Added 9_30_20 by ACL
            
            X = 1:sizex;
            Y = 1:sizey;
            figure;
            surf(X,Y,dE00,'FaceColor','interp','EdgeColor','none')
            title(['\DeltaE00 ' scanner]);
            zlim([0 100]); xlim([0 sizex]); ylim([0 sizey]);
            
        end        
        
        %% Function to run the full workflow, including:
        % acquiring images, running the sanity check, computing color truth
        % data (transmittance, XYZ, LAB, sRGB data), producing a
        % monochromatic image from the truth LAB data, and computing dE
        % between truth and mono.
        %
        function workflow_all (pathname)
            global ol490
            
            if exist(pathname,'dir') ~= 7
                mkdir(pathname);
            end
            
            path_truth = [pathname '/truth']
            path_mono = [pathname '/mono']
            path_dE = [pathname '/dE']
            
            location_filename = [pathname '/stageinfo.mat'];
            
            % MultispectralClass.find_ROI_white_black(ol490,'COM14',location_filename);
            
            if 0
                MultispectralClass.workflow_acquire(location_filename,path_truth,ol490)
            end
            
            if 0
                MultispectralClass.workflow_acquire_sanity_check(path_truth)
            end
            
            if 0
                MultispectralClass.workflow_truth(path_truth)
            end
            
            if 0
                MultispectralClass.workflow_monochrome(path_truth, path_mono)
            end
            
            if 0
                MultispectralClass.workflow_dE(path_truth, path_mono, path_dE)
            end
            
        end
        
        %% Run the workflow_all for defined sample(s)
        %
        %
        function workflow_3_tissue ()        
            %MultispectralClass.workflow_all('C:\Users\wcc\Documents\GitHub\BSC truth 11-30-2018\colon4')    
            %MultispectralClass.workflow_all('C:\Users\wcc\Documents\GitHub\BSC truth 11-30-2018\kidney4')    
            MultispectralClass.workflow_all('C:\Users\wcc\Documents\GitHub\BSC truth 11-30-2018\skin2')    
        end
        
        %% Compute dE values from HIMS data (collected from Paul's HIMS scripts)
        % Added on 9_29_20 by ACL
        % Modified version of workflow_all that skips the full workflow and
        % only computes mono and dE for the truth data (requires LAB and
        % sRGB data from the truth and sRGB data from the scanner)
        
        function workflow_dE_mono(pathname)
            global ol490
            
            if exist(pathname,'dir') ~= 7
                mkdir(pathname);
            end
            
            path_truth = [pathname '/truth']
            path_mono = [pathname '/mono']
            path_dE = [pathname '/dE']
            
            MultispectralClass.workflow_monochrome(path_truth, path_mono)
            
            MultispectralClass.workflow_dE(path_truth, path_mono, path_dE)
        end
        
        
        % Compute the dE between the scanner and the truth images
        
        
        function workflow_dE_scanner (pathname, scanner, LAB_truth, LAB_scan)
            global ol490
            
            if exist(pathname,'dir') ~= 7
                mkdir(pathname);
            end
            
            path_truth = [pathname '/truth']
            path_mono = [pathname '/mono']
            path_dE = [pathname '/dE']

            
            MultispectralClass.workflow_dE_scan(LAB_truth, LAB_scan, path_dE, scanner)
        end
        

        %% Function to show the registration and check the registration between two images
        %
        %
        function register_scanner (path_truth, path_scan, pathname, scannername)
            filepath_truth = [path_truth '/900 sRGB/truth.tif'];
            filepath_truth_lab = [path_truth '/520 CIELAB'];
            filepath_roi = [path_scan '/200 roi/scan.tif'];
            filepath_reg = [path_scan '/400 sRGB'];
            
%            MultispectralClass.register_controlpoint_save(filepath_truth,filepath_roi,filepath_reg)
%            use Matlab app to register
%            export result to movingReg
%            save movingReg to matlab_registration.mat

            [LAB_scan,LAB_truth] = MultispectralClass.register_controlpoint_show(filepath_truth,filepath_roi,filepath_reg,filepath_truth_lab);
            MultispectralClass.check_registration_with_bars(filepath_reg)
            
            % Added 9_30_20 by ACL to compute dE right after registration
            %
            % compute dE (truth to mono)
            MultispectralClass.workflow_dE_mono(pathname)
            
            % compute dE (truth to scanner)
            MultispectralClass.workflow_dE_scanner(pathname,scannername,LAB_truth,LAB_scan)

        end

        
        %% Function to define the input arguments for the register_scanner function
        %
        %
        function register_scanner_tissue
            
            datasetname = 'E:/DigitalPathology/GitHub_Repos/'
            pathname = [datasetname 'BiomaxOrgan10_Brain_H10']
            scannername = 'hamamatsu';
            
            path_truth = [pathname '/truth']
            path_scan = [pathname '/' scannername]
            
            MultispectralClass.register_scanner(path_truth,path_scan,pathname,scannername)            

        end
        
        
        %% Convert Files output from Paul's HIMS structure to MultispectralClass structure
        
        function conversion
            % Name the necessary information
            path = ('E:/DigitalPathology/GitHub_Repos/');
            path_processed = ('E:/DigitalPathology/HIMS_Data/ProcessedData/');
            Sample_Name = ('BiomaxOrgan10_UterineCervix_B10');
            Date = ('Truth');
            Scanner1 = ('hamamatsu');
            Scanner2 = ('zeiss');
            Scanner3 = ('leica');
            
            % Create directories, convert data
            fileconverter(path,path_processed,Sample_Name,Scanner1,Scanner2,Scanner3,Date)
            
        end
            
        
        % placeholder
        % does not work because of workspace
        function save_matlab_registration
            save('matlab_registration','movingReg')
        end
        
    end
    
end

