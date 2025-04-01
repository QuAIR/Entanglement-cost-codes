clear;

% Used to find irreversible (qutrit) states using EFPPT.

p = 0.0:0.001:0.025;
boundstore = zeros(4,numel(p));
for j = 1:numel(p)

    % construct states based on antisymmetric states
    rhoAB = ConstructStates(p(j));
    
    % existing bound
    Rain = rain_bound(rhoAB);
    ECW = Eeta(rhoAB);
    [LR,~] = LRbound(rhoAB);
    
    % new bound
    E_new = logfid_bineg_dual(rhoAB);

    boundstore(1, j) = ECW;
    boundstore(2, j) = LR;
    boundstore(3, j) = E_new;
    boundstore(4, j) = Rain;
end

save('boundstore_qutrit.mat', 'boundstore');


%% noisy qutrit bell
function rhoAB = ConstructStates(p)

    % Anti-symmetric state 
    v0 = [1; 0; 0];
    v1 = [0; 1; 0];
    v2 = [0; 0; 1];
    s1 = 1/sqrt(2) * (kron(v0,v1) - kron(v1,v0));
    s2 = 1/sqrt(2) * (kron(v0,v2) - kron(v2,v0));
    s3 = 1/sqrt(2) * (kron(v1,v2) - kron(v2,v1));
    rhoAB = 1*s1*s1'+1*s2*s2'+0*s3*s3';
    rhoAB = rhoAB / trace(rhoAB);
    
    rhoAB = (1-p)*rhoAB + p*eye(9)/9;
    % rhoAB = (1-p)*rhoAB + p*s3*s3';
end