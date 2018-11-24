function dE00 = lab2de00 (lab1, lab2)

    dEab = sum((lab2-lab1) .^ 2) .^ 0.5;

%-------------------------------------------------------------------
    
    Lstar1 = lab1(1);
    astar1 = lab1(2);
    bstar1 = lab1(3);
    Cstar1 = ((astar1 .^ 2) + (bstar1 .^ 2)) .^ 0.5;
    
    Lstar2 = lab2(1);
    astar2 = lab2(2);
    bstar2 = lab2(3);
    Cstar2 = ((astar2 .^ 2) + (bstar2 .^ 2)) .^ 0.5;

    dL = Lstar1 - Lstar2;
    dCab = Cstar1 - Cstar2;
    dHab = (dEab .^ 2 - dL .^ 2 - dCab .^ 2) .^ 0.5;
    Lbar = (Lstar1+Lstar2)/2;
    
    dastar = astar1 - astar2;
    dbstar = bstar1 - bstar2;

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
    SC = 1 + K1*Cstar1;
    SH = 1 + K2*Cstar1;
  
    kL = 1;
    kC = 1;
    kH = 1;

    dE94 = ((dL/kL/SL).^2 + (dCab/kC/SC).^2 + (dHab/kH/SH).^2) .^ 0.5;

%     if 0
%         kL = 1;
%         K1 = 0.045;
%         K2 = 0.015;
%     else
%         kL = 2;
%         K1 = 0.048;
%         K2 = 0.014;
%     end
% 
%     SL = 1;
%     SC = 1 + K1*Cstar1;
%     SH = 1 + K2*Cstar1;
    
    dE94 = ((dL/kL/SL).^2 + (dCab/kC/SC).^2 + (dHab/kH/SH).^2) .^ 0.5;

    
    %-------------------------------------------------------------------
    dLprime = Lstar1 - Lstar2;
    Lbar = (Lstar1+Lstar2)/2;
    Cbar = (Cstar1+Cstar2)/2;
    aprime1 = astar1 + (astar1/2) * (1 - (Cbar.^7/(Cbar.^7+25.^7)).^0.5);
    aprime2 = astar2 + (astar2/2) * (1 - (Cbar.^7/(Cbar.^7+25.^7)).^0.5);

    Cprime1 = (aprime1 .^ 2 + bstar1 .^ 2) .^ 0.5;
    Cprime2 = (aprime2 .^ 2 + bstar2 .^ 2) .^ 0.5;
    Cbarprime = (Cprime1 + Cprime2) / 2;
    dCprime = Cprime2 - Cprime1;
    
    hprime1 = mod(atan2(bstar1,aprime1)*180/pi + 360, 360);
    hprime2 = mod(atan2(bstar2,aprime2)*180/pi + 360, 360);
    
    if (abs(hprime1 - hprime2) <= 180)
        dhprime = hprime2 - hprime1;
    else if hprime2 <= hprime1
            dhprime = hprime2 - hprime1 + 360;
        else
            dhprime = hprime2 - hprime1 - 360;
        end
    end
    
    dHprime = 2 * (Cprime1 * Cprime2) .^ 0.5 * sin(dhprime*pi/180/2);
    
    if abs(hprime1 - hprime2) > 180
        Hbarprime = (hprime1 + hprime2 + 360) / 2;
    else
        Hbarprime = (hprime1 + hprime2) / 2;
    end
    
    T = 1 - 0.17*cos((Hbarprime - 30)*pi/180) + ...
        0.24 * cos(2*Hbarprime * pi/180) + ...
        0.32 * cos ((3*Hbarprime + 6 ) * pi/180) -...
        0.20 * cos((4*Hbarprime - 63) * pi/180);
    
    SL = 1 + (0.015 * (Lbar - 50) .^ 2) / ((20+(Lbar - 50).^2).^0.5);
    
    SC = 1 + 0.045 * Cbarprime;
    
    SH = 1 + 0.015 * Cbarprime * T;
    
    RT = -2 * (Cbarprime .^ 7 / (Cbarprime.^7 + 25.^7)).^0.5 * sin((60 * exp(-((Hbarprime - 275)/25).^2))*pi/180);

    dE00 = ((dLprime/kL/SL).^2 + (dCprime/kC/SC).^2 + (dHprime/kH/SH).^2 + (RT*dCprime/kC/SC*dHprime/kH/SH)) .^ 0.5;
end
