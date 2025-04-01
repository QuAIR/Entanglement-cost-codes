clear;

%% Used to find (qubit) states that EFPPT own advantages compared to 
% the existed.

% amplitude damping & depolarizing factor
pA = 0.1;
pD = 0.0:0.02:0.66;

% storing data
boundstore = zeros(5,numel(pD));
for j = 1:numel(pD)

    % qubit noisy MES state
    rhoAB = AD_Depo_qubit(pA, pD(j));

    % existing bound
    Rain = rain_bound(rhoAB);
    EoF = EntFormation(rhoAB);
    ECW = Eeta(rhoAB);
    [LR,~] = LRbound(rhoAB);
    
    % new bound
    E_new = logfid_bineg_dual(rhoAB);

    boundstore(1, j) = EoF;
    boundstore(2, j) = ECW;
    boundstore(3, j) = LR;
    boundstore(4, j) = E_new;
    boundstore(5, j) = Rain;
end

save('boundstore_qubit_fixad_varde.mat', 'boundstore');


%% Noisy bell states
function rhoAB = AD_Depo_qubit(pA, pD)

    MES = MaxEntangled(2,0,1)*MaxEntangled(2,0,1)';
    
    % Depolarizing channel
    JDepo = DepolarizingChannel(2, 1-pD);
    
    % Amplitude Damping channel
    gamma = pA;
    E1=[1 0;0 sqrt(1-gamma)];
    E2=[0 sqrt(gamma);0 0];
    JAD = ChoiMatrix({E1;E2});
    
    % Act on Bell state
    J_AD_Depo = PermuteSystems(kron(JAD,JDepo),[1 3 2 4]);
    rhoAB = ApplyMap(MES, J_AD_Depo);
end


%% Noisy bell states
function rhoAB = AD_AD_qubit(pA)

    MES = MaxEntangled(2,0,1)*MaxEntangled(2,0,1)';
    
    % Amplitude Damping channel
    gamma = pA;
    E1=[1 0;0 sqrt(1-gamma)];
    E2=[0 sqrt(gamma);0 0];
    JAD = ChoiMatrix({E1;E2});
    
    % Act on Bell state
    J_AD_AD = PermuteSystems(kron(JAD,JAD),[1 3 2 4]);
    rhoAB = ApplyMap(MES, J_AD_AD);
end