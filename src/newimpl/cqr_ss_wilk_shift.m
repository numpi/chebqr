function s = cqr_ss_wilk_shift(d,b,u,v,imax)
%CQR_SS_WILK_SHIFT 

dA = d(imax-1) * d(imax) - b(imax-1) * (conj(b(imax-1)) - conj(u(imax))*v(imax-1) + u(imax-1)*conj(v(imax)));
trA = d(imax-1) + d(imax);
D = sqrt(trA^2 - 4*dA);

l1 = trA/2 + D/2;
l2 = trA/2 - D/2;

if abs(l1 - d(imax)) < abs(l2 - d(imax))
    s = l1;
else
    s = l2;
end


end

