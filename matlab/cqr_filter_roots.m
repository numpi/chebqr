function x = cqr_filter_roots(p, x)
%CQR_FILTER_ROOTS Filter the computed roots

% We keep the root only if the inclusion radius computed using the Newton
% correction intersect the real line. 
n = length(x);
ok = zeros(1, n);

for j = 1 : length(x)
    % Compute the radius of the Bernstein Ellipse containing the given
    % root. If too large, drop it, as it is guaranteed to be outside of the
    % region of interest, and we cannot evaluate it stably anyway. 
    rho = roots([ 1, -x(j), 1 ]); 
    rho = max(abs(rho)); % Roots are coupled as g and 1/g
    
    if rho <= 2^(1 / n) && real(x(j)) >= -1 && real(x(j)) <= 1
        newt  = cqr_newton(p, x(j));
        ok(j) = abs(newt) * (n - 1) >= abs(imag(x(j)));
    end
end

x = x(ok == 1);
x = real(x);

end

