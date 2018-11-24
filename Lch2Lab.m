% convert Lch into Lab
% https://en.wikipedia.org/wiki/Lab_color_space

function lab = Lch2Lab (lch)
    lab = lch;
    if abs(lch(1,2)) < 0.000001
        lab(1,1) = lch(1,1);
        lab(1,2) = 0;
        lab(1,3) = 0;
    else
        lab(1,1) = lch(1,1);
        lab(1,2) = lch(1,2) * cos(lch(1,3)*pi/180);
        lab(1,3) = lch(1,2) * sin(lch(1,3)*pi/180);
    end
end
