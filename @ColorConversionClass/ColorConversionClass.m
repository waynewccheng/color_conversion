% Q: too many copies for the same color conversion
% A: collect all functions in one class
% 11-24-2018 Thanksgiving
% change methods to Static
classdef ColorConversionClass
    
    %ColorConversionClass used in multispectral imaging system
    %
    
    properties
        Property1
    end
    
    methods (Static)
        
        function obj = ColorConversionClass
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
        end
        
        
        
        function rgb = XYZ2sRGB (XYZ)
            %% convert a XYZ vector into sRGB
            % https://en.wikipedia.org/wiki/SRGB
            % The CIE XYZ values must be scaled so that the Y of D65 ("white")
            % is 1.0 (X,Y,Z = 0.9505, 1.0000, 1.0890). This is usually true but
            % some color spaces use 100 or other values (such as in the Lab
            % article).
            % 11-24-2018
            % 7-23-2018
            % 9-3-2015
            %        function [rgb overflow_rate underflow_rate] = XYZ2sRGB_check (XYZ)
            
            XYZ_normalized_to_one = XYZ / 100;
            
            %% constants
            m = [3.2410 -1.5374 -0.4986; -0.9692 1.8760 0.0416; 0.0556 -0.2040 1.0570];
            a = 0.055;
            
            %% linearize
            rgb = m * XYZ_normalized_to_one';
            
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
        
        %
        % CIE D65 standard illuminant
        % 11-25-2018 normalize
        %
        function ls = spd_d65 ()
            
            % to fix data path
            % spec is a 401x2 matrix
            load ('C:\Users\wcc\Documents\GitHub\color_conversion\@ColorConversionClass\datain/spec_cied65','spec');
            
            % reduce to 41x2
            ls_before_normalized = spec(1:10:401,2);
            
            % calculate CIEXYZ
            XYZ = ColorConversionClass.spd2XYZ(ls_before_normalized);
            
            % normalize CIE Y to 100
            Y = XYZ(2);
            ls = ls_before_normalized / Y * 100;
            
        end
        
        %% convert spectrum into XYZxyz
        % spd is a 41xN matrix
        % each spd is 41x1; no wavelength
        %
        function XYZ = spd2XYZ (spd)
            
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
            
            % extend x_bar from 41x1 to 41xN
            input_n = size(spd,2);
            x_bar = repmat(cmf(:,2),1,input_n);
            y_bar = repmat(cmf(:,3),1,input_n);
            z_bar = repmat(cmf(:,4),1,input_n);
            
            X = sum(spd .* x_bar);
            Y = sum(spd .* y_bar);
            Z = sum(spd .* z_bar);
            
            % output Nx3
            XYZ = [X' Y' Z'];
        end
        
        %%
        %
        % CIEXYZ to CIELAB
        %
        function ret = XYZ2lab (XYZ, XYZ_white)
            % XYZ is k-by-3
            % XYZ_white is 1-by-3
            k = size(XYZ,1);
            XYZn = repmat(XYZ_white,k,1);
            XYZ_over_XYZn = XYZ./XYZn;
            
            lstar = 116 * helpf(XYZ_over_XYZn(:,2)) - 16;
            astar = 500 * (helpf(XYZ_over_XYZn(:,1)) - helpf(XYZ_over_XYZn(:,2)));
            bstar = 200 * (helpf(XYZ_over_XYZn(:,2)) - helpf(XYZ_over_XYZn(:,3)));
            
            ret=[lstar astar bstar];
            %     if lstar > 100
            %         ['exceeding in xyz2lab']
            %         [x y z xn yn zn lstar astar bstar]
            %         lstar = 100;
            %     end
            
            return
            
            function ys = helpf (t)
                % conditional mask
                t_greater = (t > power(6/29,3));
                
                % conditional assignment
                t(t_greater) = t(t_greater) .^ (1/3);
                t(~t_greater) = t(~t_greater) * (((29/6)^2)/3) + 4/29;
                
                ys = t;
            end
        end
        
        %
        % lan1, lab2: 3xN matrix
        %
        function dEab = LAB2dEab (lab1, lab2)
            % sum in the vertical direction
            dEab = sum((lab2-lab1) .^ 2,1) .^ 0.5;
        end
        
        %
        % dE00
        % wrong results on vector ???
        % spent the whole afternoon to debug !!!
        % 11-25-2018 vectorize
        % lab1, lab2: 3xN matrix
        %
        function [dE00 dE94 dEab] = LAB2dE (lab1, lab2)
            
            dEab = sum((lab2-lab1) .^ 2,1) .^ 0.5;
            
            %-------------------------------------------------------------------
            
            Lstar1 = lab1(1,:);
            astar1 = lab1(2,:);
            bstar1 = lab1(3,:);
            Cstar1 = ((astar1 .^ 2) + (bstar1 .^ 2)) .^ 0.5;  % vector
            
            Lstar2 = lab2(1,:);
            astar2 = lab2(2,:);
            bstar2 = lab2(3,:);
            Cstar2 = ((astar2 .^ 2) + (bstar2 .^ 2)) .^ 0.5;  % vector
            
            dL = Lstar1 - Lstar2;                             % vector
            dCab = Cstar1 - Cstar2;                           % vector
            dHab = (dEab .^ 2 - dL .^ 2 - dCab .^ 2) .^ 0.5;  % vector
            Lbar = (Lstar1+Lstar2)/2;                         % vector
            
            dastar = astar1 - astar2;  % vector
            dbstar = bstar1 - bstar2;  % vector
            
            if 1
                kL = 1;
                K1 = 0.045;
                K2 = 0.015;
            else
                kL = 2;
                K1 = 0.048;
                K2 = 0.014;
            end
            
            SL = 1;
            SC = 1 + K1*Cstar1;   % vector
            SH = 1 + K2*Cstar1;   % vector
            
            kL = 1;
            kC = 1;
            kH = 1;
            
            % !!!!! vector operations !!!!!
            dE94 = ((dL./kL./SL).^2 + (dCab./kC./SC).^2 + (dHab./kH./SH).^2) .^ 0.5;
            
            
            
            %-------------------------------------------------------------------
            dLprime = Lstar1 - Lstar2;
            Lbar = (Lstar1+Lstar2)/2;
            Cbar = (Cstar1+Cstar2)/2;
            aprime1 = astar1 + (astar1/2) .* (1 - (Cbar.^7./(Cbar.^7+25.^7)).^0.5);
            aprime2 = astar2 + (astar2/2) .* (1 - (Cbar.^7./(Cbar.^7+25.^7)).^0.5);
            
            Cprime1 = (aprime1 .^ 2 + bstar1 .^ 2) .^ 0.5;
            Cprime2 = (aprime2 .^ 2 + bstar2 .^ 2) .^ 0.5;
            Cbarprime = (Cprime1 + Cprime2) / 2;
            dCprime = Cprime2 - Cprime1;
            
            hprime1 = mod(atan2(bstar1,aprime1)*180/pi + 360, 360);
            hprime2 = mod(atan2(bstar2,aprime2)*180/pi + 360, 360);
            
            %{
            if (abs(hprime1 - hprime2) <= 180)
                dhprime = hprime2 - hprime1;
            else if hprime2 <= hprime1
                    dhprime = hprime2 - hprime1 + 360;
                else
                    dhprime = hprime2 - hprime1 - 360;
                end
            end
            %}
            
            mask1 = (abs(hprime1 - hprime2) <= 180);
            mask2 = ~mask1 & (hprime2 <= hprime1);
            mask3 = ~mask1 & ~mask2;
            
            dhprime(mask1) = hprime2(mask1) - hprime1(mask1);
            dhprime(mask2) = hprime2(mask2) - hprime1(mask2) + 360;
            dhprime(mask3) = hprime2(mask3) - hprime1(mask3) - 360;
            
            dHprime = 2 * ((Cprime1 .* Cprime2) .^ 0.5) .* sin(dhprime*pi/180/2);
            
            %{
            if abs(hprime1 - hprime2) > 180
                Hbarprime = (hprime1 + hprime2 + 360) / 2;
            else
                Hbarprime = (hprime1 + hprime2) / 2;
            end
            %}
            
            mask4 = (abs(hprime1 - hprime2) > 180);
            mask5 = ~mask4;
            
            Hbarprime(mask4) = (hprime1(mask4) + hprime2(mask4) + 360) / 2;
            Hbarprime(mask5) = (hprime1(mask5) + hprime2(mask5)) / 2;
            
            T = 1 - 0.17*cos((Hbarprime - 30)*pi/180) + ...
                0.24 * cos(2*Hbarprime * pi/180) + ...
                0.32 * cos ((3*Hbarprime + 6 ) * pi/180) -...
                0.20 * cos((4*Hbarprime - 63) * pi/180);
            
            SL = 1 + (0.015 * (Lbar - 50) .^ 2) ./ ((20+(Lbar - 50).^2) .^ 0.5);
            
            SC = 1 + 0.045 * Cbarprime;
            
            SH = 1 + 0.015 * Cbarprime .* T;
            
            RT = -2 * (Cbarprime .^ 7 / (Cbarprime.^7 + 25.^7)).^0.5 * sin((60 * exp(-((Hbarprime - 275)/25).^2))*pi/180);
            
            
            dE00 = ((dLprime./kL./SL).^2 + (dCprime./kC./SC).^2 + (dHprime./kH./SH).^2 + (RT.*dCprime./kC./SC.*dHprime./kH./SH)) .^ 0.5;
            
        end
        
        %
        % Sharma's code
        % http://www2.ece.rochester.edu/~gsharma/ciede2000/dataNprograms/deltaE2000.m
        %
        function de00 = deltaE2000_Sharma (Labstd,Labsample, KLCH)
            %function de00 = deltaE2000(Labstd,Labsample, KLCH )
            % Compute the CIEDE2000 color-difference between the sample between a reference
            % with CIELab coordinates Labsample and a standard with CIELab coordinates
            % Labstd
            % The function works on multiple standard and sample vectors too
            % provided Labstd and Labsample are K x 3 matrices with samples and
            % standard specification in corresponding rows of Labstd and Labsample
            % The optional argument KLCH is a 1x3 vector containing the
            % the value of the parametric weighting factors kL, kC, and kH
            % these default to 1 if KLCH is not specified.
            
            % Based on the article:
            % "The CIEDE2000 Color-Difference Formula: Implementation Notes,
            % Supplementary Test Data, and Mathematical Observations,", G. Sharma,
            % W. Wu, E. N. Dalal, Color Research and Application, vol. 30. No. 1, pp.
            % 21-30, February 2005.
            % available at http://www.ece.rochester.edu/~/gsharma/ciede2000/
            
            de00 = [];
            
            % Error checking to ensure that sample and Std vectors are of correct sizes
            v=size(Labstd); w = size(Labsample);
            if ( v(1) ~= w(1) | v(2) ~= w(2) )
                disp('deltaE00: Standard and Sample sizes do not match');
                return
            end % if ( v(1) ~= w(1) | v(2) ~= w(2) )
            if ( v(2) ~= 3)
                disp('deltaE00: Standard and Sample Lab vectors should be Kx3  vectors');
                return
            end
            
            % Parametric factors
            if (nargin <3 )
                % Values of Parametric factors not specified use defaults
                kl = 1; kc=1; kh =1;
            else
                % Use specified Values of Parametric factors
                if ( (size(KLCH,1) ~=1) | (size(KLCH,2) ~=3))
                    disp('deltaE00: KLCH must be a 1x3  vector');
                    return;
                else
                    kl =KLCH(1); kc=KLCH(2); kh =KLCH(3);
                end
            end
            
            Lstd = Labstd(:,1)';
            astd = Labstd(:,2)';
            bstd = Labstd(:,3)';
            Cabstd = sqrt(astd.^2+bstd.^2);
            
            Lsample = Labsample(:,1)';
            asample = Labsample(:,2)';
            bsample = Labsample(:,3)';
            Cabsample = sqrt(asample.^2+bsample.^2);
            
            Cabarithmean = (Cabstd + Cabsample)/2;
            
            G = 0.5* ( 1 - sqrt( (Cabarithmean.^7)./(Cabarithmean.^7 + 25^7)));
            
            apstd = (1+G).*astd; % aprime in paper
            apsample = (1+G).*asample; % aprime in paper
            Cpsample = sqrt(apsample.^2+bsample.^2);
            Cpstd = sqrt(apstd.^2+bstd.^2);
            % Compute product of chromas and locations at which it is zero for use later
            Cpprod = (Cpsample.*Cpstd);
            zcidx = find(Cpprod == 0);
            
            
            % Ensure hue is between 0 and 2pi
            % NOTE: MATLAB already defines atan2(0,0) as zero but explicitly set it
            % just in case future definitions change
            hpstd = atan2(bstd,apstd);
            hpstd = hpstd+2*pi*(hpstd < 0);  % rollover ones that come -ve
            hpstd(find( (abs(apstd)+abs(bstd))== 0) ) = 0;
            hpsample = atan2(bsample,apsample);
            hpsample = hpsample+2*pi*(hpsample < 0);
            hpsample(find( (abs(apsample)+abs(bsample))==0) ) = 0;
            
            dL = (Lsample-Lstd);
            dC = (Cpsample-Cpstd);
            % Computation of hue difference
            dhp = (hpsample-hpstd);
            dhp = dhp - 2*pi* (dhp > pi );
            dhp = dhp + 2*pi* (dhp < (-pi) );
            % set chroma difference to zero if the product of chromas is zero
            dhp(zcidx ) = 0;
            
            % Note that the defining equations actually need
            % signed Hue and chroma differences which is different
            % from prior color difference formulae
            
            dH = 2*sqrt(Cpprod).*sin(dhp/2);
            %dH2 = 4*Cpprod.*(sin(dhp/2)).^2;
            
            % weighting functions
            Lp = (Lsample+Lstd)/2;
            Cp = (Cpstd+Cpsample)/2;
            % Average Hue Computation
            % This is equivalent to that in the paper but simpler programmatically.
            % Note average hue is computed in radians and converted to degrees only
            % where needed
            hp = (hpstd+hpsample)/2;
            % Identify positions for which abs hue diff exceeds 180 degrees
            hp = hp - ( abs(hpstd-hpsample)  > pi ) *pi;
            % rollover ones that come -ve
            hp = hp+ (hp < 0) *2*pi;
            % Check if one of the chroma values is zero, in which case set
            % mean hue to the sum which is equivalent to other value
            hp(zcidx) = hpsample(zcidx)+hpstd(zcidx);
            
            Lpm502 = (Lp-50).^2;
            Sl = 1 + 0.015*Lpm502./sqrt(20+Lpm502);
            Sc = 1+0.045*Cp;
            T = 1 - 0.17*cos(hp - pi/6 ) + 0.24*cos(2*hp) + 0.32*cos(3*hp+pi/30) ...
                -0.20*cos(4*hp-63*pi/180);
            Sh = 1 + 0.015*Cp.*T;
            delthetarad = (30*pi/180)*exp(- ( (180/pi*hp-275)/25).^2);
            Rc =  2*sqrt((Cp.^7)./(Cp.^7 + 25^7));
            RT =  - sin(2*delthetarad).*Rc;
            
            klSl = kl*Sl;
            kcSc = kc*Sc;
            khSh = kh*Sh;
            
            % The CIE 00 color difference
            de00 = sqrt( (dL./klSl).^2 + (dC./kcSc).^2 + (dH./khSh).^2 + RT.*(dC./kcSc).*(dH./khSh) );
            
            return
        end
        
        %
        %
        %
        function [dE00 dE94 dEab] = image2dE (fn1,fn2)
            
            im1 = imread(fn1);
            im2 = imread(fn2);
            
            lab1 = rgb2lab(im1,'ColorSpace','srgb','WhitePoint','d65');
            lab2 = rgb2lab(im2,'ColorSpace','srgb','WhitePoint','d65');
            
            lab1_1d = reshape(lab1,size(lab1,1)*size(lab1,2),3);
            lab2_1d = reshape(lab2,size(lab2,1)*size(lab2,2),3);
            
            [dE00_1d dE94_1d dEab_1d] = ColorConversionClass.LAB2dE(lab1_1d',lab2_1d');
            
            dE00 = reshape(dE00_1d,size(im1,1),size(im1,2));
            dE94 = reshape(dE94_1d,size(im1,1),size(im1,2));
            dEab = reshape(dEab_1d,size(im1,1),size(im1,2));
            
        end
        
        
    end
end
