function val_dual = logfid_bineg_dual(rho)
c = size(rho);
d = c(1);
cvx_begin sdp quiet
    cvx_solver sedumi
    cvx_precision best
    
    variable Q(d,d) hermitian
    variable R(d,d) hermitian
    variable U(d,d) hermitian
    variable V(d,d) hermitian
    
    t_dual = real(trace(Q*rho));
    
    minimize t_dual
    
        subject to
    
        U>=0;V>=0;
        [Q, -eye(d);
        -eye(d), R] >= 0;
        R <= PartialTranspose(U-V, 2);
        PartialTranspose(U+V,2) >= -trace(Q*rho)*eye(d);
        PartialTranspose(U+V,2) <= trace(Q*rho)*eye(d);
    
cvx_end
val_dual = -2*log2(t_dual);
end
