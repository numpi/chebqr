function [imin, imax, deflated] = cqr_ss_check_deflations(d, b, imin, imax)
%CQR_SS_CHECK_DEFLATIONS 

deflated = 0;

while imin < imax && (abs(d(imin)) + abs(d(imin+1))) * eps > abs(b(imin))
    deflated = deflated + 1;
    imin = imin + 1;
end

while imax > imin && (abs(d(imax-1)) + abs(d(imax))) * eps > abs(b(imax-1))
    deflated = deflated + 1;
    imax = imax - 1;
end

end

