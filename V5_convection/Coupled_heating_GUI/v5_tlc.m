function [rgb, cmap, Tticks] = v5_tlc(T, T_start, bw, style)
%V5_TLC  Colour a temperature field the way the thermochromic gel looks in
%        the photographs (../../50W CW 6 min; registered crops in photos/).
%
%   rgb = v5_tlc(T, T_start, bw)            T any 2-D array (degC) -> m x n x 3
%   rgb = v5_tlc(T, T_start, bw, 'vivid')   same hues, saturated (easier to read)
%   [~, cmap, Tticks] = v5_tlc([], T_start, bw, style)   256-colour bar and its temperatures
%
%   Below T_start the gel is its own hazy blue-grey.  Through the band the
%   pigment plays pink -> yellow -> green -> teal (the photographed core),
%   then blue at the top of the band, and fades back to clear above it (a
%   cholesteric pigment's clearing point).  'photo' colours are 9x9-pixel
%   means sampled from photos/photo_332s.jpg (gel, pink ring, yellow, green,
%   teal core, bluish patches), lightly de-hazed; the photographs never show
%   deep blue or a clear core, so the hottest photographed gel is still
%   inside the band.

if nargin < 4, style = 'photo'; end
s_k = [-0.06  0.00  0.12  0.35  0.55  0.78  1.00  1.12  1.35];
switch style
    case 'photo'
        gel = [0.645 0.730 0.785];
        c_k = [gel;
               0.70 0.70 0.76;             % onset: first pink specks
               0.82 0.72 0.78;             % pink ring       (sample 0.79 0.76 0.79)
               0.83 0.80 0.76;             % yellow          (0.81 0.79 0.79)
               0.70 0.86 0.78;             % green           (0.72 0.86 0.80)
               0.56 0.855 0.82;            % teal core       (0.57 0.86 0.81)
               0.50 0.72 0.86;             % blue, top of the band (not photographed)
               0.56 0.70 0.82;             % fading
               gel];                       % clear above the band
    otherwise
        gel = [0.62 0.71 0.78];
        c_k = [gel;
               0.74 0.66 0.74;
               0.88 0.62 0.70;
               0.90 0.82 0.62;
               0.66 0.86 0.64;
               0.46 0.82 0.74;
               0.30 0.58 0.86;
               0.40 0.60 0.80;
               gel];
end

if isempty(T)
    Tticks = linspace(T_start - 0.15*bw, T_start + 1.45*bw, 256);
    cmap = colour((Tticks - T_start)/bw);
    rgb = [];
    return
end
rgb = reshape(colour((T(:) - T_start)/bw), [size(T) 3]);

    function C = colour(s)
        s = min(max(s(:), s_k(1)), s_k(end));
        C = interp1(s_k, c_k, s);
    end
end
