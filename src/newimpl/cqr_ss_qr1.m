function [d,b,u,v,its,Q] = cqr_ss_qr1(d,b,u,v,aed,schurv,Q)
%CQR_SS_QR1 Single-shift QR for Hermitian plus rank 1

n = length(d);
its = 0;

imin = 1; imax = n;

last_deflation = 0;

if ~exist('aed', 'var')
    aed = true;
end

if ~exist('schurv', 'var')
    schurv = false;
end

if ~exist('Q', 'var')
    Q = zeros(n, 0);
end

ns = 0;
s = [];

ww = sort(eig(cqr2full(d,b,u,v)))

while imax > imin
    % Compute a shift for use in the next iteration
    if ns == 0
        if its <= last_deflation + 20
            s = cqr_ss_wilk_shift(d, b, u, v, imax); ns = 1;
        else
            s = rand + 1i * rand; ns = 1;
        end
    end
    
    % Compute the first bulge, and then chase it down to the end    
    [d, b, u, v, Q] = cqr_ss_chase(d, b, u, v, imin, imax, s(1), schurv, Q);
    ns = ns - 1; s = s(2:end);
    
    % Check for top / bottom deflations
    [imin, imax, deflated] = cqr_ss_check_deflations(d, b, imin, imax);
    
    if (deflated > 0)
        last_deflation = its;
    end
    
    % Check for middle deflation -- this may call the QR iteration
    % recursively if the problem is split in two parts
    if mod(its, 5) == 0
        [d, b, u, v, imin, imax, its, Q] = ...
            cqr_ss_check_middle_deflations(d, b, u, v, imin, imax, its, aed, schurv, Q);
    end
    
    % Perform aggressive early deflation, if the matrix is large enough
    if ns == 0 && aed && imax >= imin + 12
        [d, b, u, v, imin, imax, s, ns] = cqr_ss_aed(d, b, u, v, imin, imax);
    end
    
    its = its + 1;
    
    % For debugging only
    % [ sort(eig(cqr2full(d,b,u,v))), cqr2full(d,b,u,v) ]
    if n > 6
        ww = [ ww, sort(eig(cqr2full(d,b,u,v))) ]
        keyboard
    end
end

end

