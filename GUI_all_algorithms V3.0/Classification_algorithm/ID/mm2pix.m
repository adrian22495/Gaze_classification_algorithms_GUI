function xpix = mm2pix(xmm, DPI)

    MM_PER_IN = 25.4;
    xpix = xmm / (MM_PER_IN / DPI);

return