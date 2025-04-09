function val_dual = alt_logfid_bineg_dual(rho)
c = size(rho);
d = c(1);
cvx_begin sdp quiet
    cvx_solver sedumi
    cvx_precision best
    variable Q(d,d) hermitian
    variable R(d,d) hermitian
    variable S(d,d) hermitian

    loss = real(trace(rho*Q));
    minimize loss
    subject to
        [Q -eye(d);
         -eye(d) R] >= 0;
        -loss*eye(d) <= PartialTranspose(S, 2) <= loss * eye(d);
        -S <= PartialTranspose(R, 2) <= S;
cvx_end
val_dual = -2*log2(loss);
end