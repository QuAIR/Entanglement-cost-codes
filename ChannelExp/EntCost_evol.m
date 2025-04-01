clear;

dA = 2;
dB = 2^1;

dAB = dA*dB;

X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];

% type I Heisenberg model
c = [-0.5, -0.5, -1.0];
HAB = c(1) * kron(X,X) + c(2) * kron(Y,Y) + c(3) * kron(Z, Z);
JDP = full(DepolarizingChannel(dAB, 0.95));
MES = MaxEntangled(2,0,1)*MaxEntangled(2,0,1)';

% Amplitude Damping channel
gamma = 0.1;
E1=[1 0;0 sqrt(1-gamma)];
E2=[0 sqrt(gamma);0 0];
JAD = ChoiMatrix({E1;E2});
J_AD_ID = PermuteSystems(kron(JAD, dA*MES),[1 3 2 4]);

% evolution
T = 30;
T_e = 3.2;
data_store = {};

for t = 0:1/T:T_e
    Ut = expm(-1i*t*HAB);
    JUt = ChoiMatrix({Ut});
    comp = kron(eye(dAB), J_AD_ID) * kron(PartialTranspose(JUt, 2, [dAB,dAB]), eye(dAB));
    JN = PartialTrace(comp, 2, [dAB,dAB,dAB]);

    JNstate = JN / trace(JN);
    JNstateperm = PermuteSystems(JNstate, [1,3,2,4], [2,2,2,2]);

    % other bounds
    Eeta_JN = Eeta(JNstateperm);
    [Lami, ~] = LRbound(JNstateperm);
    ELami_JN = Lami;

    % our bounds
    cvx_begin sdp quiet
    % cvx_solver SDP3
    cvx_precision best
        variable Q(dAB^2, dAB^2) hermitian
        variable R(dAB^2, dAB^2) hermitian
        variable S(dAB^2, dAB^2) hermitian

        loss = trace(JNstate*Q);
        minimize loss
        subject to
            [Q -eye(dAB^2);
             -eye(dAB^2) R] >= 0;
            -loss*eye(dAB^2) <= PartialTranspose(S, [2,4], [dA, dB, dA, dB]) <= loss * eye(dAB^2);
            -S <= PartialTranspose(R, [2,4], [dA, dB, dA, dB]) <= S;
    cvx_end
    val = loss;

    data_store{end+1} = [-2*log2(val); Eeta_JN; ELami_JN];
    t
end


%% functions
% max log negativity
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
