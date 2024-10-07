function BMCC = bmcc(phi , r)
% phi: Dictionary
% r : Block size
% BMCC: Block Mutual Coherence Constant

% [1]  Y. Eldar, P. Kuppinger, and H. Bolcskei, “Block-sparse signals: Un- ¨
% certainty relations and efficient recovery,” IEEE Trans. Signal Process.,
% 2010
% [2] D. Donoho and M. Elad, “Optimally sparse representation in general
% (nonorthogonal) dictionaries via `
% 1 minimization,” in Proc. Natl. Acad.
% Sci., 2003

coeff = nan(1, 6); % Six operator norms:{1:(1,1),2:(1,2),3:(1,Inf),4:(2,2),5:(2,Inf),6:(Inf,Inf)}
for i = 1 : length(coeff) % Coefficients of the operator norm 
    coeff(i) = bmcc_coeff(r , r , r , i); % Because d is constant, it is computed outside of the following for loop
end % for i = 1 : length(coeff) % Coefficients of the operator norm 

blockNb = size(phi , 2)/r; % # of blocks
BMCC = zeros(1 , 6 + 2 + 1); % Initializing with a minimum value (six operator norms + 2 [1] + [2])
for i = 1 : blockNb % iteration over all blocks k
    blockRng = 1 : blockNb; % blocks k'
    blockRng(blockRng == i) = []; % excluding current block k
    for ii = blockRng % excluding current block; for k' != k
        % Proposed characterisation (Property III.1)
        TEMP = pinv(phi(: , (i - 1)*r + 1 : i*r)) * phi(: , (ii - 1)*r + 1 : ii*r);
        
        for k = 1 : length(coeff) % Coefficients (#=6) of the operator norm 
            temp = coeff(k) * bmcc_operatorNorm(TEMP , k);
            if temp > BMCC(k) % store the max value
                BMCC(k) = temp;
            end % if temp > BMCC(k) % store the max value
        end % for k = 1 : length(coeff) % Coefficients (#=6) of the operator norm 
 
        %% [1]'s mu_B
        temp = (1/r) * norm(phi(: , (i - 1)*r + 1 : i*r)' * phi(: , (ii - 1)*r + 1 : ii*r) , 2); % computing similarity between blocks for norm 2 (Eldar)
        if temp > BMCC(7) % store the max value
            BMCC(7) = temp;
        end % if temp > BMCC(7) % store the max value
    end % for ii = blockRng % excluding current block; for k' != k
    
    %% [1]'s sub-coherence
    temp = phi(: , (i - 1)*r + 1 : i*r)' * phi(: , (i - 1)*r + 1 : i*r);
    temp = max(abs(temp(~eye(size(temp)))));
    if temp > BMCC(8) % store the max value
        BMCC(8) = temp;
    end
end % for i = 1 : blockNb % iteration over all blocks k

%% [2]'s coherence
temp =  phi' * phi; % Gram matrix
BMCC(9) = max(abs(temp(~eye(size(temp)))));
