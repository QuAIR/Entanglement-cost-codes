clear;

%% setup

dA = 2;
dB = 2;
dAB = dA*dB;

v0 = [1;0];
v1 = [0;1];

SWAP = [1 0 0 0;
        0 0 1 0;
        0 1 0 0;
        0 0 0 1];

JN = ChoiMatrix({SWAP});
MES = 4*MaxEntangled(4, 0, 1) * MaxEntangled(4, 0, 1)';
eota = 1;

phi = (1/2)*pi;
p=0.0:0.01:1.0;
data_store = zeros(1,numel(p));

da = 2;
db = 2;
dap = 2;
dbp = 2;

for j=1:numel(p)
    
    Uphi = [1 0 0 0;
            0 exp(1i*phi) 0 0;
            0 0 exp(1i*phi) 0;
            0 0 0 exp(2i*phi)];

    K1 = sqrt(p(j)) * SWAP;
    K2 = sqrt(1-p(j)) * Uphi * SWAP;

    JNstate = kron(eye(4), K1) * MES * kron(eye(4), K1') + kron(eye(4), K2) * MES * kron(eye(4), K2');
    JNstate = JNstate / trace(JNstate);

    JNstateperm = PermuteSystems(JNstate, [1,3,2,4], [2,2,2,2]);
    data_store(2, j) = Eeta(JNstateperm);
    [Lami, ~] = LRbound(JNstateperm);
    data_store(3, j) = Lami;
    data_store(1, j) = logfid_bineg_dual(JNstateperm);
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