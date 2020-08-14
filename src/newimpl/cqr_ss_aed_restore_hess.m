function [d,b,u,v] = cqr_ss_aed_restore_hess(d,b,u,v,Q,imin,imax)
%CQR_SS_AED_RESTORE_HESS 

k = imax - imin + 1;
M = zeros(3, 3);
N = zeros(3, 3);

% Special cases first
if k == 0
    return;
end

if k == 1
    b(imax-1) = Q(1,1);
    return;
end

% The case j == k is handle separately
G = planerot(Q(k-1:k));
Q(k-1:k) = G*Q(k-1:k);

M(1,1) = d(imax-1); 
M(2,1) = b(imax-1);
M(1,2) = conj(M(2,1)) - conj(u(imax)) * v(imax-1) + u(imax-1) * conj(v(imax));
M(2,2) = d(imax);

% Update M
M(1:2,:) = G * M(1:2,:); 
M(:,1:2) = M(:,1:2) * G';

u(imax-1:imax) = G * u(imax-1:imax);
v(imax-1:imax) = G * v(imax-1:imax);

H = cqr2full(d(imin:imax), b(imin:imax-1), u(imin:imax), v(imin:imax));
H(end-1:end,end-1:end) = M(1:2,1:2);
HH{1} = H;

for j = k-1 : -1 : 2
    % Extend M
    M(3,1) = 0;
    M(3,2) = M(2,1);
    M(3,3) = M(2,2);
    M(2,1) = b(imax-k+j-1);
    M(2,2) = M(1,1);
    M(2,3) = M(1,2);
    M(1,1) = d(imax-k+j-1);
    M(1,2) = conj(M(2,1)) - conj(u(imax-k+j)) * v(imax-k+j-1) + u(imax-k+j-1) * conj(v(imax-k+j));
    M(1,3) = - conj(u(imax-k+j+1)) * v(imax-k+j-1) + u(imax-k+j-1) * conj(v(imax-k+j+1));
    
    % Determine the rotation to annihilate the entries in Q    
    G = planerot(Q(j-1:j));
    Q(j-1:j) = G * Q(j-1:j);
    
    % Apply it to M
    M(1:2,:) = G * M(1:2,:);
    M(:,1:2) = M(:,1:2) * G';    

    u(imax-k+j-1:imax-k+j) = G * u(imax-k+j-1:imax-k+j);
    v(imax-k+j-1:imax-k+j) = G * v(imax-k+j-1:imax-k+j);
    
    % Chase the bulge down to the bottom
    N = M;
    
    for i = j + 1 : k
        
    end
    
    
    % Clean the new bulge in position (3,1)
    G = planerot(M(2:3,1));
    M(2:3,:) = G * M(2:3,:);
    M(:,2:3) = M(:,2:3) * G';
        
    u(imax-k+j:imax-k+j+1) = G * u(imax-k+j:imax-k+j+1);
    v(imax-k+j:imax-k+j+1) = G * v(imax-k+j:imax-k+j+1);
    
    % Extract the computed diagonal and subdiagonal entries
    d(imax-k+j+1) = M(3,3);
    b(imax-k+j) = M(3,2);
    
    H = cqr2full(d(imin:imax), b(imin:imax-1), u(imin:imax), v(imin:imax));
    H(j-1:j+1,j-1:j+1) = M;
    sort(eig(H))
    keyboard
end

% Read the remaining entries
b(imin-1) = Q(1,1);
d(imin)   = M(1,1);
b(imin)   = M(2,1);
d(imin+1) = M(2,2);

end

