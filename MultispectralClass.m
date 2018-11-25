
classdef MultispectralClass
    
    %% Multispectral Imaging System processing
    % 11-24-2018 Thanksgiving
    
    properties
        Property1
    end
    
    methods
        
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
        function [transmittance_array, sizey, sizex] = frame2reflectance_white (obj, foldername, foldername_white)
            
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
        % D65
        % 7-23-2015
        % convert reflectance into RGB using D65
        % usage: rgb = reflectance2D65(reflectance_array7);
        % 
        function XYZ = reflectance2XYZ (obj, reflect_array, ls)
            
            disp('Combining reflectance and illuminant into XYZ...')
            
            XYZxyz0 = obj.spd2XYZ(ls);
            whiteY = XYZxyz0(2);
            
            vector_length = size(reflect_array,2);
            ls_array = repmat(ls,1,vector_length);
            spd_array = reflect_array .* ls_array;
            
            % rescale to [0:1]
            % define the white level here
            XYZ = obj.spd2XYZ(spd_array)/whiteY * 1;
            
        end
        
        %
        %
        %
        function spd_array = reflectance2spd (obj, reflect_array, ls)
            
            vector_length = size(reflect_array,2);
            ls_array = repmat(ls,1,vector_length);
            spd_array = reflect_array .* ls_array;
            
        end
        
        
        %
        %Convert the old "reflectance.mat" to "transmittance.mat"
        % 
        function convert_transmittance
            
            load('reflectance')
            sizex = 3376
            sizey = 2704
            transmittance_array = reshape(reflectance_array,41,sizey,sizex);
            
            save('transmittance','transmittance_array','sizex','sizey','-v7.3')
            
        end
        
    end
end

