% Code implementation for: 
% Fardin Afdideh, Ronald Phlypo, Christian Jutten, 
% "Exact Recovery Conditions for Optimally
% Block-sparse Representation in General Dictionaries
% via Weighted Mixed-(Pseudo-)Norm Minimization".

% SL as a function of d for different values of q and p, averaged over
% rpt=100 realisations of random dictionaries, with m=40 and n=400.

% [1]  Y. Eldar, P. Kuppinger, and H. Bolcskei, “Block-sparse signals: Un- ¨
% certainty relations and efficient recovery,” IEEE Trans. Signal Process.,
% 2010
% [2] D. Donoho and M. Elad, “Optimally sparse representation in general
% (nonorthogonal) dictionaries via `
% 1 minimization,” in Proc. Natl. Acad.
% Sci., 2003

clear
clc
close all

%% Parameters
rpt = 100; % # of repetitions
m = 40; % # of rows in the dictionary
n = 400; % # of columns in the dictionary
rNb = 6; % # of different block lengths (d)
r = find(rem(n, 1 : n) == 0, rNb, 'first'); % Chooses rNb block lengths which are divisible by the # of columns (n)

%% Load utilities (helper functions)
addpath(genpath(pwd))

%%
BMCC = nan(rpt , 6 + 2 + 1 , length(r) , 2); % Block Mutual Coherence Constant (BMCC): repetitions x methods x d x regular/ort
rLen = length(r);
rng('default') % The same random numbers are produced as if you restarted MATLAB.
t0 = tic;
parfor i = 1 : rpt % Iterate over rpt repetitions
    phiIni = randn(m, n);
    for j = 1 : rLen % block size
        phi = phi_generator(phiIni, n/r(j));
        for k = 1 : 2 % length(phi): regular/ort
            BMCC(i , : , j , k) = bmcc(phi{k} , r(j));
            %         disp(['',num2str(i),'/',num2str(rpt),' - ',num2str(j),'/',num2str(length(r)),' - ',num2str(k),'/',num2str(length(phi)),''])
        end % for k = 1 : 2 % regular/ort
    end % for j = 1 : length(r) % block size
end % for i = 1 : rpt % Iterate over rpt repetitions
elapsedTime = toc(t0);
clear t0 i j k phiIni phi

%% Block-sparse exact recovery conditions (BERC) computation
BERC_mean = nan(size(BMCC , 3) , size(BMCC , 2)-1 , size(BMCC , 4));
BERC_std = nan(size(BERC_mean));
for i = 1 : size(BMCC , 3) % Different block size
    for j = 1 :  size(BMCC , 4) % Regular - orthonormal dictionary
        % Proposed BERC (Theorem III.3)
        BERC = r(i).*((1 + (r(i).*BMCC(: , 1:6 , i , j)).^(-1)) ./ 2);
        BERC_mean(i , 1:6 , j) = mean(BERC); % across rpt
        BERC_std(i , 1:6 , j) = std(BERC); % across rpt
        
        % [1]'s BERC
        mu_B = BMCC(: , 7 , i , j);
        nu = BMCC(: , 8 , i , j);
        BERC = r(i).*((1 + (r(i).*mu_B).^(-1) - (((r(i) - 1) .* nu ) ./ (r(i) .* mu_B) ) ) ./ 2);
        BERC_mean(i , 7 , j) = mean(BERC); % across rpt
        BERC_std(i , 7 , j) = std(BERC); % across rpt
        
        % [2]'s ERC
        ERC = (1 + BMCC(: , 9 , i , j).^(-1) ) ./ 2;
        BERC_mean(i , 8 , j) = mean(ERC); % across rpt
        BERC_std(i , 8 , j) = std(ERC); % across rpt
    end
end
clear i j BERC ERC mu_B nu

%% Fig. 1
rSel = 1:size(BERC_mean, 1); % select desired block size (r) indices to plot
tmpMean = BERC_mean(rSel , : , 1);
tmpMeanOrth = BERC_mean(rSel , : , 2);
tmpStd = BERC_std(rSel , : , 1);
tmpStdOrth = BERC_std(rSel , : , 2);
tmpStd(tmpMean < 0) = 0;
tmpMean(tmpMean < 0) = 0;
tmpStdOrth(tmpMeanOrth < 0) = 0;
tmpMeanOrth(tmpMeanOrth < 0) = 0;

figure
% ort dictionary
h1 = barweb(tmpMeanOrth, tmpStdOrth , [], r, [], 'Block size (r)', 'Sparsity level (SL)', [], [], [], [], []);
drawnow
grayColor = ones(1, 3); % White bars (ort)
for i = 1 : length(h1.bars)
    h1.bars(i).FaceColor = grayColor;
    h1.bars(i).LineWidth = 1;
end

hold on
% regular dictionary
h2 = barweb(tmpMean, tmpStd, [], r, [], 'Block size (r)', 'Sparsity level (SL)', [], [], [], [], []);
extraCol = 2; % edge colors to be discarded
grayColor = gray(size(tmpMean, 2) + 2*extraCol);
grayColor = grayColor(extraCol+1:size(grayColor, 1)-extraCol, :);
for i = 1 : length(h2.bars)
    h2.bars(i).FaceColor = grayColor(i, :);
    h2.bars(i).LineWidth = 1;
end

set(gca,'XTick', 1 : size(tmpMean, 1))
legendTxt = {'M_{1,1}({\bf\Phi})' , 'M_{1,2}({\bf\Phi})' , 'M_{1,\infty}({\bf\Phi})' , ...
    'M_{2,2}({\bf\Phi})' , 'M_{2,\infty}({\bf\Phi})' ,...
    'M_{\infty,\infty}({\bf\Phi})' , '[5]' , '[10]', '{\bf\Phi}_{ort}'};
curveSel = get(gca, 'Children');
legend(curveSel([length(curveSel)-2*size(BERC_mean, 2):-1:length(curveSel)-3*size(BERC_mean, 2)+1, length(curveSel)]), legendTxt, 'location', 'west')
set(gca,'YTick', 1 : 12, 'YLim', [0 11])
grid on
axis square

clear dSel tmpMean tmpMeanOrth tmpStd tmpStdOrth groupnames bw_xlabel h1 h2 i extraCol grayColor legendTxt curveSel