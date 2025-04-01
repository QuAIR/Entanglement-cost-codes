%% Eeta from Wang et. al.
function out = Eeta(rhoAB)
    d = max(size(rhoAB));
    da = sqrt(d);
    db = da;
    [V, D] = eig(rhoAB);
    P = 0;
    for i=1:d
        if D(i,i) > 0.00001
            P = P+V(:,i)*V(:,i)';
        end
    end
    cvx_begin sdp quiet
        variable Y(d,d) hermitian
        
        Yb = PartialTranspose(Y, 2, [da db]);
        f = norm(Yb, inf);
    
        minimize f
        subject to
            PartialTranspose(P, 2, [da db]) + Y >= 0;
            PartialTranspose(P, 2, [da db]) - Y <= 0;
    cvx_end
    out = -log2(f);
end

