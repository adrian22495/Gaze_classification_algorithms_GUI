function [ medianBCEA ] = medianBCEA68( test )
%Input: Información del test obtenida mediante store_stim_info o load_and_compile
%outpu: Media entre las fBCEA de todos los estimulos del test
fBCEA = test(8:20:size(test, 2));
fBCEA(fBCEA == -1) = []
medianBCEA = median(fBCEA);
end

