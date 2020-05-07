function vf = filterNoiseSGolay(data, validity, k, f)
% k and f are the Salvitzky-Golay filtering parameters
if nargin < 4
    k = 3;
    f = 5;
end

vf = [];
Xforf = [];

if (validity(1))
    Xforf = data(1);
else
    vf = data(1);
end

for i=2:size(data,1)
    if (validity(i))   % in valid span
        Xforf = [Xforf; data(i)];
    else        % in gap
        if (validity(i-1)) % start of gap
            % Filter
            if (size(Xforf,1) >= f)
                XforfFilt = sgolayfilt(Xforf,k,f);
            else
                XforfFilt = Xforf;
            end
            % Copy to output
            vf = [vf; XforfFilt];
            % Set Xforf to empty for next valid span
            Xforf = [];
        end
        % I copy the zeros (if it is the start of gap, this is only done
        % after filtering and copying the valid span)
        vf(i,1) = data(i);
    end
end
% Since I only copy at the start of gap, I need to copy the possible last
% valid span here
if (size(Xforf,1) >= f)
    XforfFilt = sgolayfilt(Xforf,k,f);
else
    XforfFilt = Xforf;
end
% Copy to output
vf = [vf; XforfFilt];

return
