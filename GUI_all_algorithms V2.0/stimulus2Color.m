function c = stimulus2Color(id, n, V)
    if nargin < 3
        V = 1;
    end
    S = 1;
    H = 0.67*(1 - id/(n-1));
    c = hsv2rgb([H S V]);
end