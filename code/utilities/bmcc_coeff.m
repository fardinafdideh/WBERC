function coeff = bmcc_coeff(r_1 , r_2 , r_max , qpCounter)
% r_1: Length of first block (r_k in the manuscript)
% r_2: Length of second block (r_k' in the manuscript)
% r_max: Maximum length of blocks
% qpCounter: {1:(1,1),2:(1,2),3:(1,Inf),4:(2,2),5:(2,Inf),6:(Inf,Inf)}
% coeff: The coefficient in the formulum of BMCC (Property III.1)

if qpCounter == 1  % q=1; p=1
    coeff = r_2 / (r_max * r_1);
elseif qpCounter == 2 % q=1; p=2 
    coeff = r_2 / (r_max * (r_1^(1/2)));
elseif qpCounter == 3 % q=1; p=inf
    coeff = r_2 / r_max;
elseif qpCounter == 4 % q=2; p=2
    coeff = (r_2^(1/2)) / (r_max * (r_1^(1/2)));
elseif qpCounter == 5 % q=2; p=inf
    coeff = (r_2^(1/2)) / r_max;
elseif qpCounter == 6 % q=inf; p=inf
    coeff = 1 / r_max;
end
