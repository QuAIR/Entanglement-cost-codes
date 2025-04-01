function val = ED_werner_holevo_ub(d)
% compute the upper bound of distillable entanglement of Werner-Holevo
% channel
if mod(d,2) == 0;
    val = log2((d+2)/d);
else 
    val = 0.5*log2((d+3)/(d-1));
end
val = val;
end

