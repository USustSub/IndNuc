function [wm,w0,wp]=centred3(hm,hp)
%CENTRED3 Three-point centred first-derivative weights on unequal spacing.
%
%   [wm,w0,wp] = CENTRED3(hm,hp) gives
%
%       f'(s) ~ wm*f(s-hm) + w0*f(s) + wp*f(s+hp)
%
%   exact for quadratics and therefore second order for ANY hm, hp. With
%   hm == hp == h it returns (-1, 0, +1)/(2h), the plain centred difference,
%   so the centre weight vanishes and the operator is unchanged on a uniform
%   grid.
%
%   WHY THE CENTRE WEIGHT MATTERS. The plain (f+ - f-)/(hm+hp) form has
%
%       (f+ - f-)/(hm+hp) = f' + (hp-hm)/2 * f'' + O(h^2)
%
%   and for a smooth map hp-hm = dzeta^2*s'', so it is still second-order
%   CONVERGENT but its leading error is proportional to s'', the curvature of
%   the mesh grading. These weights are exact for quadratics, so that term is
%   identically absent.
%
%   Now a wrapper over FDWEIGHTS.

w=fdweights(0,[-hm 0 hp],1);
wm=w(1); w0=w(2); wp=w(3);
end
