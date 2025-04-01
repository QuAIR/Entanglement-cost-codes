clear;

%%
sample = 500;
d = 3;
r = 2;
rho_i = RandomDensityMatrix(d^2, 0, r);

bound_store = zeros(4, sample);
for j=1:sample
    UR = RandomUnitary(d^2);
    rho_f = UR * rho_i * UR';
    E_new = logfid_bineg_dual(rho_f);
    Rains = rain_bound(rho_f);
    E_eta = Eeta(rho_f);
    [LR, ~] = LRbound(rho_f);

    bound_store(1, j) = Rains;
    bound_store(2, j) = E_new;
    bound_store(3, j) = E_eta;
    bound_store(4, j) = LR;

    j
end 

