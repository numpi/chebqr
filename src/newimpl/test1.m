n = 16;
rng(0);
d = randn(n, 1);
b = randn(n-1, 1);
u = rand(n,1);
v = rand(n,1);

if n <= 4096
    tic;
    H = cqr2full(d,b,u,v); 
    ee = eig(H, 'nobalance'); ee = sort(ee);
    toc
    %E = randn(n); E = E / norm(E, 'fro');
    %ee2 = eig(H + norm(H, 'fro') * eps * E); ee2 = sort(ee2);
    ee2 = ee;
else
    ee = zeros(n,1);
    ee2 = ee;
end

if n <= 512
    tic;
    [d1, b1, u1, v1, its1, Q1] = cqr_ss_qr1(d, b, u, v, true, true, eye(n));
    toc
    its1
else
    d1 = zeros(n,1);
end

tic; 
[d2, its2] = cqr_ss_qr1_mex(d, b, u, v, true);
toc
its2

return;

sort([ ee, d1, d2 ])
close all
plot(real(ee), imag(ee), 'bo')
hold on; plot(real(ee2), imag(ee2), 'ko');
hold on; plot(real(d1), imag(d1), 'rx')
hold on; plot(real(d2), imag(d2), 'g+')