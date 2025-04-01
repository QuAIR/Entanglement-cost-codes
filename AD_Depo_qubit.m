%----------------------------------------------------------------------------------
%                    The Bell state destroyed by the
%     Qubit Amplitude Damping Channel and Qubit Depolarizing Channel
%----------------------------------------------------------------------------------
% Required packages: 
% CVX http://cvxr.com/cvx/download/
% QETLAB http://www.qetlab.com/Main_Page
%----------------------------------------------------------------------------------
% Standard usage: rhoAB = AD_Depo_qubit(pA, pD)
% Variables:
%        pA     -     noise parameter of the amplitude damping channel
%                       e.g. pA = 0.1
%        pD     -     noise parameter of the depolarizing channel
%                       e.g. pD = 0.22
%----------------------------------------------------------------------------------

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
J_AD_Depo = PermuteSystems(kron(JDepo, JAD),[1 3 2 4]);
rhoAB = ApplyMap(MES, J_AD_Depo);
end
