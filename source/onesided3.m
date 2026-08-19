function [w1,w2,w3]=onesided3(h1,h2)
%ONESIDED3 Three-point one-sided first-derivative weights at the first node.
%
%   [w1,w2,w3] = ONESIDED3(h1,h2) gives
%
%       f'(s1) ~ w1*f(s1) + w2*f(s1+h1) + w3*f(s1+h1+h2)
%
%   second order for ANY h1, h2. With h1 == h2 == h it returns
%   (-3, 4, -1)/(2h), the uniform weights the free-surface rows in build_LH
%   carried before the mesh could be stretched -- which is why the uniform
%   regression against build_LH_reference still holds.
%
%   Now a wrapper over FDWEIGHTS, which is the single source of truth for
%   every non-uniform stencil in the operator. Kept as a named function
%   because the free-surface rows read better with it.

w=fdweights(0,[0 h1 h1+h2],1);
w1=w(1); w2=w(2); w3=w(3);
end
