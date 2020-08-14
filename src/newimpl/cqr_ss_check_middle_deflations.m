function [d, b, u, v, imin, imax, its, Q] = cqr_ss_check_middle_deflations(d, b, u, v, imin, imax, its, aed, schurv, Q)
%CQR_SS_CHECK_MIDDLE_DEFLATIONS 

di = abs(d(imin));

for i = imin + 1 : imax - 2
     % This avoids computing absolute values twice
     dim1 = di;
     di = abs(d(i));
     
     if (abs(b(i)) < (di + dim1) * eps)
        [d(imin:i), b(imin:i-1), u(imin:i), v(imin:i), lits, Q(imin:i,:)] = ...
            cqr_ss_qr1(d(imin:i), b(imin:i-1), u(imin:i), v(imin:i), aed, schurv, Q(imin:i,:));
        its = its + lits;
        imin = i + 1;
     end
end

end

