function ksi=build_ksi(G,L,dy,a,b,sigman0)
%BUILD_KSI Lapusta-style adaptive-timestep coefficient, one per fault node.
%   dy may be a scalar or a per-node vector of the LOCAL down-dip spacing.
%   It must be the local value: on a y-stretched mesh part of the
%   rate-and-state fault can sit in the coarsened arm, and feeding the core
%   spacing everywhere understates the cell size there. That errs safe (too
%   small dy -> larger k1 -> smaller ksi -> smaller dt) but wastes steps and
%   is not the stiffness the arm actually has.
    k1=pi/4*G./dy.*L./a./sigman0;
    k2=(b-a)./a;
    k3=(k1-k2).^2/4-k1;
    k4=min(1./(k1-k2), 0.2);
    k5=min(1-k2./k1, 0.2);
    ksi=k4.*(k3>0)+k5.*(~(k3>0));
end