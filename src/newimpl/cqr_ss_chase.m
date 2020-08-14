function [d, b, u, v, Q] = cqr_ss_chase(d, b, u, v, imin, imax, s, schurv, Q)
%CQR_SS_CHASE

if (imin == imax)
    return;
end

% Find the first rotation to use
G = planerot([ d(imin) - s; b(imin) ]);

% Construct the relevant block, and apply the rotation
M = [ d(imin) , conj(b(imin)) - conj(u(imin+1)) * v(imin) + u(imin) * conj(v(imin+1)) ; ...
      b(imin) , d(imin+1) ; ...
      0    , 0 ];
  
if imax > imin + 1
    M(3,2) = b(imin+1);
end
  
M(1:2,:) = G * M(1:2,:);
M = M * G';

u(imin:imin+1) = G * u(imin:imin+1);
v(imin:imin+1) = G * v(imin:imin+1);

Q = cqr_update_schur_vectors(G, imin, schurv, Q);

d(imin) = M(1,1);

for i = imin : imax - 2
     % Find the rotation required to clean the bulge
     G = planerot(M(2:3,1));
     
     % Update the first column of M, and store the updated
     % subdiagonal entry; note that we only update from the left, as the
     % right transformation will leave this column unchanged. 
     w = G * M(2:3,1); b(i) = w(1);
     
     % Create the new M, and apply G
     M(1,1) = M(2,2);
     M(2,1) = M(3,2);
     M(1,2) = conj(M(2,1)) - conj(u(i+2)) * v(i+1) + u(i+1) * conj(v(i+2));
     M(2,2) = d(i+2);
     
     if i < imax - 2
        M(3,2) = b(i+2);
        M(3,1) = 0;
     end

     M(1:2,:) = G * M(1:2,:);
     M = M * G';
     
     d(i+1) = M(1,1);

     u(i+1:i+2) = G * u(i+1:i+2);
     v(i+1:i+2) = G * v(i+1:i+2);
     
     Q = cqr_update_schur_vectors(G, i+1, schurv, Q);
end

% At the end, we need to copy the trailing entries
d(imax) = M(2,2);
b(imax-1) = M(2,1);

end

