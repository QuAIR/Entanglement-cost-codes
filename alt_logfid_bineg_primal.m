function val_primal = alt_logfid_bineg_primal(rho)
c = size(rho);
d = c(1);
cvx_begin sdp quiet
    cvx_solver sedumi
    cvx_precision best
    variable X(d,d) complex
    variable M(d,d) complex semidefinite
    variable N(d,d) complex semidefinite
    variable C(d,d) complex semidefinite
    variable D(d,d) complex semidefinite

    loss = 0.5*real(trace(X + X'));
    maximize loss
    subject to
        [rho, X;
          X', PartialTranspose(C - D,2)] >= 0;
        PartialTranspose(M - N, 2) == C + D;
        trace(M+N) <= 1;
cvx_end
val_primal = -2*log2(loss)
end