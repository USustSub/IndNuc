function w=fdweights(x0,p,m)
%FDWEIGHTS Finite-difference weights on arbitrarily spaced nodes.
%
%   w = FDWEIGHTS(x0,p,m) returns row vector w with
%
%       f^(m)(x0)  ~  sum_i w(i) * f(p(i))
%
%   exact for every polynomial up to degree numel(p)-1, for ANY node positions
%   p and ANY evaluation point x0. x0 need NOT be one of the nodes -- half-node
%   stencils are the reason this function exists.
%
%   Method: impose exactness on the monomials, i.e. solve A*w = e_m with
%   A(k+1,i) = (p(i)-x0)^k / k!. Bengt Fornberg (Math. Comp. 51 (1988) 699-706)
%   gives a stable O(n^2) recursion for the same weights; at the n = 3 and 4
%   used here the direct solve is the same answer in three lines. Nodes are
%   scaled to unit spread before the solve, which keeps the Vandermonde-like
%   system well conditioned when p spans kilometres.
%
%   This is the single source of truth for every derivative stencil on a
%   non-uniform mesh. ONESIDED3 and CENTRED3 are thin wrappers over it and
%   return exactly the textbook uniform weights when the spacings are equal:
%
%       fdweights(0,[0 h 2h],1)*2h   ->  (-3,  4, -1)      = onesided3
%       fdweights(0,[-h 0 h],1)*2h   ->  (-1,  0,  1)      = centred3
%       fdweights(0,[-1.5 -.5 .5 1.5]*h,2)*2h^2 -> (1,-1,-1,1)
%
%   The last one is the half-node second derivative that build_LH previously
%   hard-coded. On the nys = 36 seam spacing (100, 100, 116.5 m) the correct
%   weights are NOT antisymmetric -- the two middle weights differ by more than
%   a factor of two -- and the hard-coded set is 25 % wrong on a pure quadratic.

p=p(:).';
n=numel(p);
if m>=n
    error('fdweights:tooFewNodes', ...
        'need at least %d nodes for a derivative of order %d; got %d.', ...
        m+1,m,n);
end
d=p-x0;
s=max(abs(d));
if s==0
    error('fdweights:degenerate','all nodes coincide with the evaluation point.');
end
q=d/s;                       % scaled to unit spread: conditioning
A=zeros(n);
for k=0:n-1
    A(k+1,:)=q.^k/factorial(k);
end
b=zeros(n,1);
b(m+1)=1;
w=(A\b).'/s^m;               % undo the scaling: d/dx = (1/s) d/dq
end
