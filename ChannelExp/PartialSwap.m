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
eota = 1i;

p=0.0:0.01:1.0;
data_store = zeros(1,numel(p));

da = 2;
db = 2;
dap = 2;
dbp = 2;

for j=1:numel(p)
    
    Up = sqrt(p(j))*eye(4) + eota*sqrt(1-p(j))*SWAP;
    JNstate = ChoiMatrix({Up});
    JNstate = JNstate / trace(JNstate);

    JNstateperm = PermuteSystems(JNstate, [1,3,2,4], [2,2,2,2]);
    data_store(1, j) = logfid_bineg_dual(JNstateperm);

    data_store(2, j) = Eeta(JNstateperm);
    [Lami, ~] = LRbound(JNstateperm);
    data_store(3, j) = Lami;
    j

end
