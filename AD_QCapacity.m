function val = AD_QCapacity(eta)
% compute the quantum capacity of qubit amplitude damping channel
delta = 0.00001;
Q = [];
if eta <= 0.5
    val = 0;
elseif 0.5 < eta <= 1
    for p = 0.01:delta:1.0
        Q(end+1) = -eta*p*log2(eta*p) - (1-eta*p)*log2(1-eta*p) + ...
            (1-eta)*p*log2((1-eta)*p) + (1-(1-eta)*p)*log2((1-(1-eta)*p));
    end 
    val = max(Q);
else
    val = 0;
    error('Wrong eta');
end
val = val;
end

