function Q = cqr_update_schur_vectors(G, i, schurv, Q)
%CQR_UPDATE_SCHUR_VECTORS

if schurv
    Q(i:i+1,:) = G * Q(i:i+1,:);
end

end

