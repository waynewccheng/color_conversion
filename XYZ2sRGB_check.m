%% convert a XYZ vector into sRGB
% 11-24-2018 Thanksgiving
% 7-23-2018
% 9-3-2015

function [rgb overflow_rate underflow_rate] = XYZ2sRGB_check (XYZ)

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

% 7-25-2018
% highlight overflow- and underflow-pixels
if 0
    
    overflow_rate = nnz(rgb > 1) / size(rgb,2)
    underflow_rate = nnz(rgb < 0) / size(rgb,2)
    
    %% check overflow
    overflow_mask = (rgb(1,:) > 1) | (rgb(2,:) > 1) | (rgb(3,:) > 1);
    rgb(1,overflow_mask) = 1;
    rgb(2,overflow_mask) = 0;
    rgb(3,overflow_mask) = 0;
    
    %% check uderflow
    underflow_mask = (rgb(1,:) < 0) | (rgb(2,:) < 0) | (rgb(3,:) < 0);
    rgb(1,underflow_mask) = 0;
    rgb(2,underflow_mask) = 1;
    rgb(3,underflow_mask) = 0;
    
else
    
    overflow_rate = 0;
    underflow_rate = 0;
    
end

%% comply with the old form
rgb = rgb';

end
