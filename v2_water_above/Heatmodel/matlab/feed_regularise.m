function [q, info] = feed_regularise(q, xg, yg, zg, Tc)
%FEED_REGULARISE  Optionally remove the feed-point singularity from the
%                 microwave power deposition.
%
%   [q, info] = FEED_REGULARISE(q, xg, yg, zg, Tc)
%
%   WHY THERE IS A SINGULARITY
%   The field model prescribes the antenna current I(x) = I_m sin(k_e(L-x))
%   on 0 <= x <= L and zero outside.  With L = lambda_e/4 that current is
%   MAXIMUM at the feed, so truncating it there is a step discontinuity in I,
%   and continuity of charge, dI/dx + j w rho = 0, puts a point charge
%       Q = I_m / w = 4.02 / 1.539e10 = 2.6e-10 C
%   at (x,y,z) = (0,0,h).  Its quasi-static field 0.5 mm away in a medium with
%   eps' = 74.3 is Q/(4 pi eps0 eps' r^2) = 1.3e5 V/m, which is exactly the
%   magnitude the solver returns there.  At the TIP the current goes smoothly
%   to zero, so there is no such charge and no such spike -- which is the
%   check that this is the mechanism.
%
%   In a real applicator the feed current does not stop: it continues onto the
%   outside of the coaxial outer conductor.  The field model does not include
%   that sheath (its Sec. 11, limitation 3), so the deposition within a couple
%   of millimetres of the feed is not trustworthy.  A real feed junction IS a
%   hot spot -- the "backward heating along the shaft" reported in the
%   catheter-ablation literature -- but neither its magnitude nor its shape is
%   predicted correctly here.
%
%   THE BRAID ARM REPLACES ALL OF THIS.  MODEL_CONFIG's L_sheath now carries
%   the return current on the outside of the coax braid, which is what the
%   current at the feed actually does -- it returns, it does not stop and it
%   does not pile up as charge.  With a continuous I(x) there is no delta in
%   dI/dx, no point charge, and nothing to remove, so GUI_MODEL sets
%   Tc.r_feed_excl = 0 and this function becomes the identity: no mask, no grey
%   semicircle in the maps, and the reported peak is a real prediction instead
%   of the value on the edge of an arbitrary sphere.
%
%   Everything below therefore applies only when L_sheath = 0, i.e. to the
%   v1/v2 comparison case.
%
%   WHAT THIS FUNCTION DOES
%     Tc.feed_fix = 'none'  (default, and the BASELINE)
%         q is returned unchanged.  The thermal model is then exactly the
%         field model plus conduction, with no hidden edits, and it agrees
%         with E_gel everywhere.  The reported PEAK temperature excludes a
%         sphere of radius Tc.r_feed_excl about the feed, because inside it
%         the answer is meaningless and outside it the answer is conservative.
%     Tc.feed_fix = 'cap'
%         q inside that sphere is capped at the largest value found on the
%         shell just outside it.  This removes the point charge, makes the
%         peak temperature converge under grid refinement, and moves it to a
%         sensible place under the radiating element.  It also discards the
%         ~15 W the point charge was depositing, which is why it is not the
%         baseline: that power is not fictitious, we simply do not know where
%         it really goes.
%
%   Away from the feed the two agree to better than 1 %: at (6, 0, -2) mm the
%   steady-state rise is 939.9 K raw and 932.4 K capped.

% r_feed_excl = 0 means "there is no singularity to hide" -- the braid arm has
% already removed it.  Use a STRICT inequality so that the sphere is then
% genuinely empty; with <= it would still swallow the single node sitting
% exactly on the feed, which is enough to punch a visible hole in the map and
% to drop the true peak out of the reported maximum.
R = sqrt((reshape(xg,[],1,1) - Tc.feed_xyz(1)).^2 + ...
         (reshape(yg,1,[],1) - Tc.feed_xyz(2)).^2 + ...
         (reshape(zg,1,1,[]) - Tc.feed_xyz(3)).^2);
inside = R < Tc.r_feed_excl;
shell  = R >= Tc.r_feed_excl & R < Tc.r_feed_excl + 4*Tc.dh;

info = struct('mask_peak', ~inside, 'n_inside', nnz(inside));
info.P_inside_before = local_power(q, inside, Tc);

switch lower(Tc.feed_fix)
    case 'none'
        info.qcap = NaN;
    case 'cap'
        if ~any(shell(:))
            error('feed_regularise:shell','the exclusion shell is empty');
        end
        info.qcap = max(q(shell));
        q(inside) = min(q(inside), info.qcap);
    otherwise
        error('feed_regularise:mode','unknown feed_fix ''%s''', Tc.feed_fix);
end
info.P_inside_after = local_power(q, inside, Tc);
end

function P = local_power(q, m, Tc)
% trapezoidal control-volume weights are 1 in the interior; the feed sphere
% never touches a boundary face except the surface, where the weight is 1/2
w = ones(size(q));  w(:,:,1) = 0.5;  w(:,:,end) = 0.5;
P = sum(w(m).*q(m))*Tc.dh^3;
end
