function ksi=build_ksi(G,L,dy,a,b,sigman0)
    k1=pi/4*G/dy.*L./a./sigman0;
    k2=(b-a)./a;
    k3=(k1-k2).^2/4-k1;
    k4=min(1./(k1-k2), 0.2);
    k5=min(1-k2./k1, 0.2);
    ksi=k4.*(k3>0)+k5.*(~(k3>0));
end