function mm = deg2mm(degrees, viewDistance)
if nargin < 2
    viewDistance = 500;
end
mm = 2 * viewDistance * tan(deg2rad(degrees)/2);
