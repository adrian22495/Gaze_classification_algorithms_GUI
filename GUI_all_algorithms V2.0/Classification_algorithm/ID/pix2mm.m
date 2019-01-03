function xm = pix2mm(xpix, DPI)

    MM_PER_IN = 25.4;
    xm = xpix * MM_PER_IN / DPI;

return