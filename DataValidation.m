%% Evaluation R_ij and R_0 for Covid-19 in Hunan province via DBM and NGM
clear; clc; close all;

omega1 = 0.1429;
omega2 = 0.1429;
gamma1 = 0.2;
gamma2 = 0.2;
f = 3.552e-3;
p = 0.2015;
kappa = 1;


% matrices of beta_ij
Beta_med = [1.70e-13   1.12e-12    4.10e-10    1.66e-11;...
            7.16e-12   1.24e-9     1.21e-16    1.45e-10;...
            6.63e-11   1.46e-11    1.64e-19    7.71e-9; ...
            7.20e-14   3.61e-10    3.07e-9     3.66e-17    ]; % the entry-wise median of Beta matrix
 
Beta_1 = [3.62e-16  1.43e-11    5.15e-13    3.69e-14;...
          3.06e-13  8.97E-09    2.58e-15    1.41e-11;...
          8.38e-12  1.51e-17    5.34e-20    5.43e-9; ...
          3.46E-11	5.88E-11	9.81E-08    2.83E-17       ];
  


% population size of Hunan and Jilin Province
N_hunan = [13618898   26623844    20035661    8709900];
N_jilin = [3291955, 10827458, 9525131, 3395414];
P = N_hunan / sum(N_hunan);

Beta_med = Beta_med.*N_hunan;
Beta_1 = Beta_1 .* N_hunan;


c = kappa*omega2*p/(gamma2*(omega1 - omega1*p + omega2*p)) + omega1*(1 - p)/((f + gamma1)*(omega1 - omega1*p + omega2*p));

% matrix of R_ij
Rmataver = Beta_med*c
Rmat1 = Beta_1*c



R0_NGM_1 = max(eig(Beta_1)) * c
R0_DEF_1 = sum(P'.*sum(Rmat1,2))


%% For all time segments from 1 to 4
% read Betas
betaTable = readtable('Betas.xlsx');
betaTable(:,1) = [];
names = betaTable.Properties.VariableNames;

% population
N_hunan = [13618898   26623844    20035661    8709900];
P_hunan = N_hunan / sum(N_hunan);

% add max beta
beta_max = max(betaTable{:,:});
betaTable = [betaTable; mat2cell(beta_max,1,ones(1,16))];


%%% compute R_ij matrix and R0 on Beta matrices for time segments from 1 to 4
% from time segment 1 to 4
h1 = figure;
for i = 1:4
    % extract Beta matrix B from table
    B = zeros(4);
    for k = 1:16
        ii = eval(names{k}(2));
        jj = eval(names{k}(3));
        B(ii,jj) = betaTable{i,k}; 
    end
    
    %%% compute the R_ij matrix and R0 using DBM and NGM
    % compute the constant coefficient c (which is identical in DBM and NGM, due to the adoption of the group-irrelevant parameters )
    c = kappa*omega2*p/(gamma2*(omega1 - omega1*p + omega2*p)) + omega1*(1 - p)/((f + gamma1)*(omega1 - omega1*p + omega2*p));  

    % R_ij matrix for time segment i
    Rmat = B .* N_hunan * c;
    
    % R0 for segment i 
    R0_DEF = sum(P_hunan'.*sum(Rmat,2));    % definition-based method
    R0_NGM = max(eig(B.*N_hunan)) * c;      % next-generation method
    

    %%% visualizing the R_ij matrix by heatmaps 
    subplot(2,2,i);     % plot the heatmap for time segment i 
    format shorte;
    xvalue = {'0 to 14','15 to 44','45 to 64', '\geq 65'};
    yvalue = {'0 to 14','15 to 44','45 to 64', '\geq 65'};
    h_DEF = heatmap(xvalue,yvalue,Rmat);
    cellDate = {'Jan. 5 2021', 'Jan. 25 2021', 'Jan. 31 2021', 'Feb. 5 2021', 'Feb. 19 2021'};

%   h_DEF.Title = ['R_{ij} matrix, time segment ',num2str(i),' from ',cellDate{i},' to ',cellDate{i+1}];
    h_DEF.Title = ['R_{ij} Matrix of Time Segment ',num2str(i)];
    h_DEF.XLabel = 'to age group j';        % to the entirely susceptiable age group j
    h_DEF.YLabel = 'from age group i';      % from one infected individual in age group i

    h_DEF.FontName = 'Times New Roman';
    h_DEF.CellLabelFormat = '%.4e ';
%     c = colorbar;
%     c.Label.String = 'interactive R_{ij} between different groups';
    
    % save record
    R0_DEFs(i,1) = R0_DEF;
    R0_NGMs(i,1) = R0_NGM;
    eval(['B',num2str(i),'=','B;']);
    eval(['Rmat',num2str(i),'=','Rmat;']);
    eval(['R0_DEF',num2str(i),'=','R0_DEF;']);
    eval(['R0_NGM',num2str(i),'=','R0_NGM;']);
end

% heat map of increasement of Rij matrices (from time segment 1 to 4)
h2 = figure;
for i = 2:4    
    subplot(3,1,i-1);
    format shorte;
    xvalue = {'0 to 14','15 to 44','45 to 64', '\geq 65'};
    yvalue = {'0 to 14','15 to 44','45 to 64', '\geq 65'};
    h_DEF = heatmap(xvalue,yvalue,eval([' Rmat',num2str(i),' - Rmat',num2str(i-1)]));
    h_DEF.Title = ['The matrix of increasement of R_{ij}, time segment ',num2str(i-1),' to sement ',num2str(i)];
    h_DEF.XLabel = 'In the entirely susceptiable age group j';
    h_DEF.YLabel = 'One infected individual in age group i';
    h_DEF.CellLabelFormat = '%.4e ';
end


% print(h1,'-dpdf','h1.pdf','-r300');
% print(h2,'-dpdf','h2.pdf','-r300');


