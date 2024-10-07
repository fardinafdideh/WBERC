function operatorNorm = bmcc_operatorNorm(mat , qpCounter)
% mat : A matrix that its qpCounter operator norm is going to be computed
% qpCounter: {1:(1,1),2:(1,2),3:(1,Inf),4:(2,2),5:(2,Inf),6:(Inf,Inf)}
% operatorNorm: Operator Norm of the matrix "mat"

if qpCounter == 1  % q=1; p=1
    operatorNorm = norm(mat , 1);
elseif qpCounter == 2 % q=1; p=2
    operatorNorm = max(sqrt(sum(abs(mat).^2,1))); % Maximum ?2 norm of a column
elseif qpCounter == 3 % q=1; p=inf
    operatorNorm = max(max(abs(mat))); % Maximum ?? norm of a column
elseif qpCounter == 4 % q=2; p=2
    operatorNorm = norm(mat , 2);
elseif qpCounter == 5 % q=2; p=inf
    operatorNorm = max(sqrt(sum(mat.^2,2))); % Maximum ?2 norm of a row
elseif qpCounter == 6 % q=inf; p=inf
    operatorNorm = norm(mat , 'inf');
end
