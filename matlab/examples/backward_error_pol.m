function e = backward_error_pol(p, x, d)
%BACKWARD_ERROR_POL Backward error for the numerical root X of P
%
% E = BACKWARD_ERROR_POL(P, X) computes an upper to the backward error for
% the computed roots X(j) of the polynomial defined by its coefficients in 
% the Chebyshev basis:
%
%  P_1 T_0(X(i)) + P_2 T_1(X(i)) + ... + P_{n+1} T_n(X(i)),  i = 1, ..., n
%
% The backward error is defined as the minimum of ||P - Q||_2, where Q is a
% scalar multiple of the vector of Chebyshev coefficients of 
%
%  \prod_{i=1}^n (x - X(i))

n = length(x);
y = cos( (n:-1:0) / n * pi).';

p = p(:);
x = x(:);
y = y(:);

% Evaluate the polynomial at the Chebsyhev points
if d <= 16
    w = ones(n+1, 1);
    for j = 1 : n
        w = w .* (y - x(j));
    end
    
else
    bits = 256;
    w = backward_error_pol_mp(p(end:-1:1), x, bits);
end

% Make sure w is a column vector
w = w(:);

% Interpolate Q
v = [ w(end:-1:1); w(2:n) ];

%if d <= 16
    q = real(fft(v)) / n;
%else
%    ru = exp(1i * (0 : 2*n-1) ./ n * pi);
%    V = vander(ru);
%    q = real(V' * v / n); q = q(end:-1:1);
%end

q = [q(1) / 2; q(2:n); q(n+1) / 2 ];
q = q / norm(q);
p = p / norm(p);

if length(q) < length(p)
    q(length(p)) = 0;
end

% Find the projection of q along p
r = dot(p, q);

% Compute the residual
e = norm(p - r * q);

end