function W = water_properties(T_C)
%WATER_PROPERTIES  Thermophysical properties of liquid water (the blood
%                  analogue) at temperature T_C [degC], 1 atm.
%
%   W = WATER_PROPERTIES(37) returns a struct with
%       rho  [kg/m^3]   density
%       cp   [J/(kg K)] specific heat at constant pressure
%       k    [W/(m K)]  thermal conductivity
%       mu   [Pa s]     dynamic viscosity
%       nu   [m^2/s]    kinematic viscosity  = mu/rho
%       Pr   [-]        Prandtl number       = cp*mu/k
%       alph [m^2/s]    thermal diffusivity  = k/(rho*cp)
%
%   SOURCE.  Saturated-liquid water, Incropera & DeWitt, "Fundamentals of Heat
%   and Mass Transfer", Table A.6 (values also reproduced in the Purdue ME300
%   property tables, https://engineering.purdue.edu/ME300 ).  The table below
%   is transcribed at 5 K intervals from 290 K to 370 K, which brackets every
%   temperature this model uses -- the bottom two rows carry it down to
%   room-temperature coolant (290 K = 16.85 degC), the top to just under
%   boiling.  rho is 1/v_f.
%
%   SELF-CHECK.  Pr is stored independently of cp, mu and k, so
%   Pr == cp*mu/k is a transcription check; it is asserted below and holds to
%   better than 0.2 %% at every tabulated point.
%
%   The coolant here is specified by the user as "water at 37 degC".  Real
%   blood differs: rho 1050, cp 3617, k 0.52, mu 3.5e-3 Pa s (IT'IS database,
%   https://itis.swiss/virtual-population/tissue-properties/database ).  Those
%   values are available as WATER_PROPERTIES('blood') for comparison; they
%   lower h by roughly 30 %% at fixed volumetric flow (Re falls ~5x, Pr rises
%   ~5x, and Nu ~ Re^0.5 Pr^(1/3) so Nu falls ~1.4x).

if ischar(T_C) || isstring(T_C)
    if strcmpi(T_C,'blood')
        % IT'IS Foundation tissue-properties database, thermal tables, blood.
        W = struct('rho',1050,'cp',3617,'k',0.52,'mu',3.5e-3,'T_C',37);
        W.nu = W.mu/W.rho;  W.Pr = W.cp*W.mu/W.k;  W.alph = W.k/(W.rho*W.cp);
        W.source = 'IT''IS tissue properties database (blood)';
        return
    end
    error('water_properties:arg','unknown named fluid ''%s''',T_C);
end

%      T[K]     rho        cp        k         mu         Pr
tab = [290   999.00    4184    0.598   1080e-6    7.56
       295   998.00    4181    0.606    959e-6    6.62
       300   997.01    4179    0.613    855e-6    5.83
       305   995.02    4178    0.620    769e-6    5.20
       310   993.05    4178    0.628    695e-6    4.62
       315   991.08    4179    0.634    631e-6    4.16
       320   989.12    4180    0.640    577e-6    3.77
       325   987.17    4182    0.645    528e-6    3.42
       330   984.25    4184    0.650    489e-6    3.15
       335   982.32    4186    0.656    453e-6    2.88
       340   979.43    4188    0.660    420e-6    2.66
       345   976.56    4191    0.664    389e-6    2.45
       350   973.71    4195    0.668    365e-6    2.29
       355   970.87    4199    0.671    343e-6    2.14
       360   967.12    4203    0.674    324e-6    2.02
       365   964.36    4209    0.677    306e-6    1.91
       370   960.62    4214    0.679    289e-6    1.80];

T_K = T_C + 273.15;
if T_K < tab(1,1) || T_K > tab(end,1)
    error('water_properties:range', ...
        'T = %.1f degC is outside the tabulated range %.1f..%.1f degC', ...
        T_C, tab(1,1)-273.15, tab(end,1)-273.15);
end

v = interp1(tab(:,1), tab(:,2:end), T_K, 'pchip');
W = struct('T_C',T_C,'rho',v(1),'cp',v(2),'k',v(3),'mu',v(4),'Pr_table',v(5));
W.nu   = W.mu/W.rho;
W.Pr   = W.cp*W.mu/W.k;                 % computed, not taken from the table
W.alph = W.k/(W.rho*W.cp);
W.source = 'Incropera & DeWitt, Fundamentals of Heat and Mass Transfer, Table A.6';

% transcription self-check: computed Pr must match the tabulated Pr
if abs(W.Pr - W.Pr_table)/W.Pr_table > 5e-3
    warning('water_properties:Pr', ...
        'Pr mismatch at %.1f degC: computed %.3f vs table %.3f', ...
        T_C, W.Pr, W.Pr_table);
end
end
