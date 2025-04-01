%% Rains bound
function rain = rain_bound(rho)
    c = size(rho);
    d = c(1);
    cvx_begin sdp quiet
    variable K(d,d) hermitian;
    t1 = (quantum_rel_entr(rho,K)/log(2));
    minimize t1
        subject to
            K >= 0;
            SchattenNorm(PartialTranspose(K,2),1) <= 1;
    cvx_end
    rain = t1;
end

