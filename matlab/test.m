% This is just a simple test to check if everything works

n = 256;
rng(0);
d = rand(n,1); beta=rand(n-1,1);
u = rand(n,1); v = rand(n,1);

tic;
x = cqr_qr_ss_aed_par(d, beta, u, v);
toc;

tic;
y = cqr_qr_ss_aed(d, beta, u, v);
toc;

tic;
A = diag(d) + diag(beta,-1);
H = A - u*v'; H = tril(H) + tril(H,-1)';
A = H + u*v';
z = eig(A);
toc;

[ sort(x), sort(y), sort(z) ]
[ sort(abs(x)) - sort(abs(z)), sort(abs(y)) - sort(abs(z)) ]