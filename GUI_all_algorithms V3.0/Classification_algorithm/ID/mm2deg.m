function xdeg = mm2deg(xmm,vd)

    xdeg = 2 * atan(xmm / 2 / vd) * 360 / 2 / pi;
    
return
    %X = atan(X/test.VD) * 360 / PI; % -CHECK if using atan is ok in the general case
    %Y = atan(Y/test.VD) * 360 / PI; % -CHECK if using atan is ok in the general case