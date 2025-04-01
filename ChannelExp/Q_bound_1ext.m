% This script computes the upper bound in Eq . (8) of arXiv:1912.00931 .
% Required packages: 
% CVX http://cvxr . com/cvx/download/
% QETLAB http://www . qetlab . com/Main_Page
% CVX_quad https://github . com/hfawzi/cvxquad

%% input the depolarizing channel
function q_bound = Q_bound_1ext(rho)
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
        -loss*eye(d) <= PartialTranspose(S, 2, [sqrt(d), sqrt(d)]) <= loss * eye(d);
        -S <= PartialTranspose(R, 2, [sqrt(d), sqrt(d)]) <= S;
cvx_end
q_bound=loss;
end

