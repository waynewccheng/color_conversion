function theta = xy2deg (x, y)
    if abs(x) < 0.000001
        if abs(y) < 0.000001
            theta = 0;
        elseif y > 0 
            theta = 90;
        else
            theta = 270;
        end
    else
        theta = atan(abs(y)/abs(x))*180/pi;
        if x >= 0 && y >= 0 
            theta = theta;
        elseif x < 0 && y < 0
            theta = theta + 180;
        elseif x >= 0 && y < 0 
            theta = 360 - theta;
        else
            theta = 180 - theta;
        end
    end
end
