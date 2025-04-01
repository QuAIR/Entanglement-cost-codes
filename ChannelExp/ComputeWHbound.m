clear;

%%
d_range = [5];
EDBound = zeros(1,numel(d_range));
ECBound = zeros(1,numel(d_range));

for j = 1:numel(d_range)
    JW = WHChoi(d_range(j));
    EDBound(j) = ED_werner_holevo_ub(d_range(j));
    rho = JW;
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


%% Choi state of Werner-Holevo channel
function JW = WHChoi(d)
    mes = MaxEntangled(d)*MaxEntangled(d)';
    % JW = (1/(d*(d-1))) * (eye(d^2) -  SWAPOP(d));
    JW = (1/(d*(d-1))) * (eye(d^2) -  PartialTranspose(mes,2)*d);
end

function F = SWAPOP(d)
    data = eye(d);
    F = 0;
    for i = 1:d
        for j = 1:d
            F = F + kron(data(:, i) * data(:, j)', data(:, j) * data(:, i)');
        end
    end
    F = F;
end

