%--------------------------------------------------------------------------
%          Generate Choi Matrix of the Random Mixed Unitary Channel
%--------------------------------------------------------------------------
% Required packages: 
% CVX http://cvxr.com/cvx/download/
% QETLAB http://www.qetlab.com/Main_Page
%--------------------------------------------------------------------------
% Standard usage: JU = RandMixedU(p)
% Variables:
%        p     -     noise parameters of the mixed unitary channel
%                    e.g. p = [0.1, 0.2, 0.3, 0.4]
%--------------------------------------------------------------------------


function JU = RandMixedU(p)

d = 3;

p1=p(1);
p2=p(2);
p3=p(3);
p4=p(4);

U1 = RandomUnitary(d);
U2 = RandomUnitary(d);
U3 = RandomUnitary(d);
U4 = RandomUnitary(d);

E1=sqrt(p1)*U1;
E2=sqrt(p2)*U2;
E3=sqrt(p3)*U3;
E4=sqrt(p4)*U4;

JU = ChoiMatrix({E1;E2;E3;E4});

end