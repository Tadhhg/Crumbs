function [key, fields] = v5_key(P)
%V5_KEY  Exact key of the fields of P that change the stored responses (v5_solve).
%        Two P with the same key share a solution; the other fields only
%        rescale it (v5_eval).  The key is the hex of each value's bits, so it
%        is exact and needs no Java (was an MD5 via java.security, which fails
%        on MATLAB installs without a Java runtime).
fields = {'f','L','a','b','eps_ins','l_f','d', 'eps_g1','eps_g2','eps_w1','eps_w2', ...
          'rho','c','k', 'h_s','h_w','L_c', 'gel_L','gel_H','gel_W','x_start', ...
          'dim','dh','t_end','dt_frame'};
v = cellfun(@(f) P.(f), fields);
key = sprintf('%bx', v);
end
