classdef ColorConversionClass
    
    %ColorConversionClass used in multispectral imaging system
    %
    
    properties
        Property1
    end
    
    methods
        
        function obj = ColorConversionClass
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
        end
        
       
        %% convert a XYZ vector into sRGB
        % 11-24-2018
        % 7-23-2018
        % 9-3-2015
        %        function [rgb overflow_rate underflow_rate] = XYZ2sRGB_check (XYZ)
        
        function rgb = XYZ2sRGB (obj, XYZ)
            
            %% constants
            m = [3.2410 -1.5374 -0.4986; -0.9692 1.8760 0.0416; 0.0556 -0.2040 1.0570];
            a = 0.055;
            
            %% linearize
            rgb = m*XYZ';
            
            %% conditional mask
            rgb_lessorequal = (rgb <= 0.0031308);
            
            %% conditional assignment
            rgb(rgb_lessorequal) = rgb(rgb_lessorequal) * 12.92;
            rgb(~rgb_lessorequal) = (1+a)*(rgb(~rgb_lessorequal).^(1/2.4)) - a;
            
            %             % 7-25-2018
            %             % highlight overflow- and underflow-pixels
            %             if 0
            %
            %                 overflow_rate = nnz(rgb > 1) / size(rgb,2)
            %                 underflow_rate = nnz(rgb < 0) / size(rgb,2)
            %
            %                 %% check overflow
            %                 overflow_mask = (rgb(1,:) > 1) | (rgb(2,:) > 1) | (rgb(3,:) > 1);
            %                 rgb(1,overflow_mask) = 1;
            %                 rgb(2,overflow_mask) = 0;
            %                 rgb(3,overflow_mask) = 0;
            %
            %                 %% check uderflow
            %                 underflow_mask = (rgb(1,:) < 0) | (rgb(2,:) < 0) | (rgb(3,:) < 0);
            %                 rgb(1,underflow_mask) = 0;
            %                 rgb(2,underflow_mask) = 1;
            %                 rgb(3,underflow_mask) = 0;
            %
            %             else
            %
            %                 overflow_rate = 0;
            %                 underflow_rate = 0;
            %
            %             end
            
            %% comply with the old form
            rgb = rgb';
            
        end
        
                
        function ls = spd_d65 (obj)
            % to fix data path
            load ('C:\Users\wcc\Documents\GitHub\color_conversion\datain/spec_cied65','spec');
            ls = spec(1:10:401,2);
        end
        
        % convert spectrum into XYZxyz
        function XYZ = spd2XYZ (obj, s)
            
            % CIEXYZ 1931
            cmf = [
                380.0 0.001368 0.000039 0.006450;
                390.0 0.004243 0.000120 0.020050;
                400.0 0.014310 0.000396 0.067850;
                410.0 0.043510 0.001210 0.207400;
                420.0 0.134380 0.004000 0.645600;
                430.0 0.283900 0.011600 1.385600;
                440.0 0.348280 0.023000 1.747060;
                450.0 0.336200 0.038000 1.772110;
                460.0 0.290800 0.060000 1.669200;
                470.0 0.195360 0.090980 1.287640;
                480.0 0.095640 0.139020 0.812950;
                490.0 0.032010 0.208020 0.465180;
                500.0 0.004900 0.323000 0.272000;
                510.0 0.009300 0.503000 0.158200;
                520.0 0.063270 0.710000 0.078250;
                530.0 0.165500 0.862000 0.042160;
                540.0 0.290400 0.954000 0.020300;
                550.0 0.433450 0.994950 0.008750;
                560.0 0.594500 0.995000 0.003900;
                570.0 0.762100 0.952000 0.002100;
                580.0 0.916300 0.870000 0.001650;
                590.0 1.026300 0.757000 0.001100;
                600.0 1.062200 0.631000 0.000800;
                610.0 1.002600 0.503000 0.000340;
                620.0 0.854450 0.381000 0.000190;
                630.0 0.642400 0.265000 0.000050;
                640.0 0.447900 0.175000 0.000020;
                650.0 0.283500 0.107000 0.000000;
                660.0 0.164900 0.061000 0.000000;
                670.0 0.087400 0.032000 0.000000;
                680.0 0.046770 0.017000 0.000000;
                690.0 0.022700 0.008210 0.000000;
                700.0 0.011359 0.004102 0.000000;
                710.0 0.005790 0.002091 0.000000;
                720.0 0.002899 0.001047 0.000000;
                730.0 0.001440 0.000520 0.000000;
                740.0 0.000690 0.000249 0.000000;
                750.0 0.000332 0.000120 0.000000;
                760.0 0.000166 0.000060 0.000000;
                770.0 0.000083 0.000030 0.000000;
                780.0 0.000042 0.000015 0.000000;
                ];
            
            % show color matching functions
            if 0
                clf
                hold on;
                plot(cmf(:,1),cmf(:,2),'Color','r');
                plot(cmf(:,1),cmf(:,3),'Color','g');
                plot(cmf(:,1),cmf(:,4),'Color','b');
                title('Color Matching Functions')
            end
            
            input_n = size(s,2);
            x_bar = repmat(cmf(:,2),1,input_n);
            y_bar = repmat(cmf(:,3),1,input_n);
            z_bar = repmat(cmf(:,4),1,input_n);
            
            %    k = 100/10.6858;
            k = 1;
            
            X = k * sum(s .* x_bar);
            Y = k * sum(s .* y_bar);
            Z = k * sum(s .* z_bar);
            
            XYZ = [X' Y' Z'];
        end
        
    end
end

