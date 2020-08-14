function H = cqr2full(d,b,u,v,bulge)
%CQR2FULL 

if ~exist('bulge', 'var')
    bulge = 0;
end

A = diag(d) + diag(b,-1) - u*v';
H = diag(d) + diag(b,-1) + triu(A', 1) + triu(u*v',1);

end

