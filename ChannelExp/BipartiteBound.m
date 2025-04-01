%% setup
clear;

dA = 2;
dB = 2;
dAB = dA*dB;

% set-up ABA'B'
% rho = RandomDensityMatrix(dAB);
CX = [1 0 0 0;
      0 1 0 0;
      0 0 0 1;
      0 0 1 0];

SWAP = [1 0 0 0;
        0 0 1 0;
        0 1 0 0;
        0 0 0 1];

JN = ChoiMatrix({SWAP});
MES = dA*MaxEntangled(dA, 0, 1) * MaxEntangled(dA, 0, 1)';

p=0.0:0.01:1.0;
data_store = zeros(2,numel(p));

for j=1:numel(p)
    
    % Depolarizing channel
    JDepAB = DepolarizingChannel(dAB,1-p(j))*eye(dAB^2);
    Jnoise = JDepAB;

    JNTB = PartialTranspose(JN, 2, [dAB,dAB]);
    JNcomp = PartialTrace(kron(JNTB, eye(dAB))*kron(eye(dAB), Jnoise), 2, [dAB,dAB,dAB]);
    JNstate = JNcomp/trace(JNcomp);

    cvx_begin sdp quiet
    % cvx_solver SDP3
    cvx_precision best
        variable Q(dAB^2,dAB^2) hermitian
        variable R(dAB^2,dAB^2) hermitian
        variable S(dAB^2,dAB^2) hermitian

        loss = trace(JNstate*Q);
        minimize loss
        subject to
            [Q -eye(dAB^2);
             -eye(dAB^2) R] >= 0;
            -loss*eye(dAB^2) <= PartialTranspose(S, [2,4], [dA, dB, dA, dB]) <= loss * eye(dAB^2);
            -S <= PartialTranspose(R, [2,4], [dA, dB, dA, dB]) <= S;
    cvx_end
    val = loss;

    data_store(1, j) = -2*log2(val);
    j

end

%% max log negativity
function val = MaxLogNeg(JN, dim)

    dA = dim(1);
    dB = dim(2);
    dAB = dA*dB;

    % l1
    cvx_begin sdp quiet
    cvx_precision best
    cvx_solver sedumi

    variable P(dAB^2, dAB^2) hermitian %AiBiAoBo
    variable k
    
    P_TB = PartialTranspose(P, [2 4], [dA, dB, dA, dB]);
    PAB = PartialTrace(P, [3 4], [dA, dB, dA, dB]);

    minimize k
    subject to
        -P_TB <= PartialTranspose(JN, [2 4], [dA, dB, dA, dB]) <= P_TB;
        P >= 0;
        -k*eye(dAB) <= PAB <= k*eye(dAB);
    cvx_end

    l1 = log2(k);

    % l2
    cvx_begin sdp quiet
    cvx_solver sedumi
    
        variable P(dAB^2, dAB^2) hermitian %AiBiAoBo
        variable k
        
        P_TB = PartialTranspose(P, [2 4], [dA, dB, dA, dB]);
        PAB = PartialTrace(P, [3 4], [dA, dB, dA, dB]);
        minimize k
        subject to
            -P_TB <= PartialTranspose(JN, [2 4], [dA, dB, dA, dB]) <= P_TB;
            -k*eye(dAB) <= PartialTranspose(PAB, 2, [dA, dB]) <= k*eye(dAB);
            P >= 0;
    cvx_end
    
    l2 = log2(k);

    val = max([l1, l2]);
end 