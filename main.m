clc;
close all;
% 

R = 5; % Number of independent simulations, for nuclear receptors (4) 
       % we recommond to set the number of simulation up to 20 or more 
       % since the performance is not so stable because the network 
       % is too small
K = 10; % K-fold cross validation

% Initialize parameter alpha 
% Adjust these parameters will effect the performances of the method 
% The optimal parameters for AUC may differ from those for AUPR
paraT = [0.13 0.12 0.18 0.26 0.1]; % alpha for target similarity A  
paraD = [0.23 0.15 0.15 0.25 0.15];  % alpha for drug similarity A^T
paraST= [1.0 1.2 1.4 1.5]; % alpha for target similarity matrix
paraSD= [1.0 1.0 1.0 1.0]; % alpha for drug similarity matrix

data = 1 % Dataset to run, change this value for running different datasets

if data == 1
    dg = dlmread('dataset/en_sim_dg.txt');
    dc = dlmread('dataset/en_sim_dc.txt');
    adj = dlmread('dataset/en_adj.txt');
elseif data == 2
    dg = dlmread('dataset/ic_sim_dg.txt');
    dc = dlmread('dataset/ic_sim_dc.txt');
    adj = dlmread('dataset/ic_adj.txt');
elseif data == 3
    dg = dlmread('dataset/gpcr_sim_dg.txt');
    dc = dlmread('dataset/gpcr_sim_dc.txt');
    adj = dlmread('dataset/gpcr_adj.txt');
elseif data == 4
    dg = dlmread('dataset/nr_sim_dg.txt');
    dc = dlmread('dataset/nr_sim_dc.txt');
    adj = dlmread('dataset/nr_adj.txt');
elseif data == 5
    net = dlmread('dataset/matador_adj.txt');
    adj = construct_DTI_matrix(net);
end

precision_SD = []; recall_SD = [];
precision_ST = []; recall_ST = [];
precision_AT = []; recall_AT = [];
precision_AD = []; recall_AD = [];
precision_A = []; recall_A = [];
precision_ADT = []; recall_ADT =[];

aucSD = zeros(1,K);
auprSD = zeros(1,K);
aucST = zeros(1,K);
auprST = zeros(1,K);
aucAT = zeros(1,K);
auprAT = zeros(1,K);
aucAD = zeros(1,K);
auprAD = zeros(1,K);
aucA = zeros(1,K);
auprA = zeros(1,K);
aucADT = zeros(1,K);
auprADT = zeros(1,K);

auc_SD = zeros(1,R);
aupr_SD = zeros(1,R);
auc_ST = zeros(1,R);
aupr_ST = zeros(1,R);
auc_AT = zeros(1,R);
aupr_AT = zeros(1,R);
auc_AD = zeros(1,R);
aupr_AD = zeros(1,R);
auc_A = zeros(1,R);
aupr_A = zeros(1,R);
auc_ADT = zeros(1,R);
aupr_ADT = zeros(1,R);

for r = 1 : R % Number of simulations 
    disp('===========================================================');
    disp(['============= r = ' num2str(r) ' ========== of ' num2str(R) ' =======']);
    disp('===========================================================');
    fprintf('\n');
    
    y = adj;
    %[tr te] = divideMatrixCV(adj,K); % Dividing the dataset to k-fold subsets  
    %crossval_idx = divideMatrixCR(adj,K,2); % Dividing the dataset to k-fold subsets 
    crossval_idx = crossvalind('Kfold',y(:),K);
        
    for i = 1 : K % each fold at a time
        disp(['--- Run ' num2str(r) ' of ' num2str(R) ', k ' num2str(i) ' of ' num2str(K) ' ---']);

        
        train_idx = find(crossval_idx~=i);
        test_idx  = find(crossval_idx==i);

        y_train = y;
        y_train(test_idx) = 0;
        train = y_train;
        test = y - train;
        
        yy=y;
        yy(yy==0)=-1;

        if data == 5 % Working with Matador 
        
            % Note that for Matador, the target similarity obtained from 
            % the interaction information plays more important role than 
            % the drug similarity and the even the combination 

            X_AT = solve_lrr(train',train', paraT(data));  % X^*_AT
            Z_AT = (train'*X_AT)';
            statsST = evaluate_performance(Z_AT(test_idx),yy(test_idx),'classification');
            aucST(i) = statsST.auc;
            auprST(i) = statsST.aupr;
            
            X_AD = solve_lrr(train,train, paraD(data));  % X^*_A
            Z_AD = train*X_AD;
            statsSD = evaluate_performance(Z_AD(test_idx),yy(test_idx),'classification');
            aucSD(i) = statsSD.auc;
            auprSD(i) = statsSD.aupr;
            
            Z_A = (Z_AT + Z_AD)/2;
            statsA = evaluate_performance(Z_A(test_idx),yy(test_idx),'classification');
            aucA(i) = statsA.auc;
            auprA(i) = statsA.aupr;
        else
      
            % Working with drug compound 
            dc = (dc+dc')/2;
            X_SD = solve_lrr(dc,dc, paraSD(data)); % Computing X^*_{SD}
            Z_D = train * X_SD;      % Projecting A onto  X^*_{SD}, Z_D
                 
            statsSD = evaluate_performance(Z_D(test_idx),yy(test_idx),'classification');
            aucSD(i) = statsSD.auc;
            auprSD(i) = statsSD.aupr;
            
            % Working with protein sequence (target)
            dg = (dg + dg')/2;  
            X_ST = solve_lrr(dg,dg, paraST(data));  % Compute X^*_{ST}
            Z_T = X_ST * train;      % Projecting A onto X^*_{ST}, Z_T
            statsST = evaluate_performance(Z_T(test_idx),yy(test_idx),'classification');
            aucST(i) = statsST.auc;
            auprST(i) = statsST.aupr;
                 
            % Working with adjacency matrix (interaction information),
            % target similarity 
            X_AT = solve_lrr(train,train,paraT(data));  % Compute X^*_{AT}
            Z_AT = train * X_AT;   
            statsAT = evaluate_performance(Z_AT(test_idx),yy(test_idx),'classification');
            aucAT(i) = statsAT.auc;
            auprAT(i) = statsAT.aupr;
            
            % Working with adjacency matrix (interaction information),
            % drug similarity 
            X_AD = solve_lrr(train',train',paraD(data));       % Compute X^*_{AD} 
            Z_AD = (train' * X_AD)';                  % Projecting A on X^*_{D}
            statsAD = evaluate_performance(Z_AD(test_idx),yy(test_idx),'classification');
            aucAD(i) = statsAD.auc;
            auprAD(i) = statsAD.aupr;

    
            % Combine the two Z_AT and Z_AD
            Z_A = (Z_AT + Z_AD)/2;
            statsA = evaluate_performance(Z_A(test_idx),yy(test_idx),'classification');
            aucA(i) = statsA.auc;
            auprA(i) = statsA.aupr;             

            Z_ADT = 0.25*Z_T + 0.25*Z_D + 0.5*Z_A; % Combine all the information
            
            statsADT = evaluate_performance(Z_ADT(test_idx),yy(test_idx),'classification');
            aucADT(i) = statsADT.auc;
            auprADT(i) = statsADT.aupr;
        
        end
        
        fprintf('\n');
    end
    fprintf('\n');
    if data ~= 5
        auc_SD(r) = mean(aucSD);
        aupr_SD(r) = mean(auprSD);
    
        auc_ST(r) = mean(aucST);
        aupr_ST(r) = mean(auprST);

        auc_AT(r) = mean(aucAT);
        aupr_AT(r) = mean(auprAT);
        
        auc_AD(r) = mean(aucAD);
        aupr_AD(r) = mean(auprAD);

        auc_A(r) = mean(aucA);
        aupr_A(r) = mean(auprA);
    
        auc_ADT(r) = mean(aucADT);
        aupr_ADT(r) = mean(auprADT);
        
    else
        auc_AT(r) = mean(aucAT);
        aupr_AT(r) = mean(auprAT);
      
        auc_AD(r) = mean(aucAD);
        aupr_AD(r) = mean(auprAD);
        
        auc_A(r) = mean(aucA);
        aupr_A(r) = mean(auprA);
    end

end

%% Print results 
if data == 5
    disp(['The AUC from Z_ST: ' num2str(auPRAT(auc_AT))]);
    disp(['The AUPR from Z_ST: ' num2str(mean(aupr_AT))]);
    
    disp(['The AUC from Z_SD: ' num2str(auPRAT(auc_AD))]);
    disp(['The AUPR from Z_SD: ' num2str(mean(aupr_AD))]);
    
    disp(['The AUC from Z_A: ' num2str(auPRAT(auc_A))]);
    disp(['The AUPR from Z_A: ' num2str(mean(aupr_A))]);
else
    disp(['The AUC from Z_SD: ' num2str(mean(auc_SD))]);
    disp(['The AUPR from Z_SD: ' num2str(mean(aupr_SD))]);
    fprintf('\n');
    
    disp(['The AUC from Z_ST: ' num2str(mean(auc_ST))]);
    disp(['The AUPR from Z_ST: ' num2str(mean(aupr_ST))]);
    fprintf('\n');
    
    disp(['The AUC from Z_AD: ' num2str(mean(auc_AT))]);
    disp(['The AUPR from Z_AD: ' num2str(mean(aupr_AT))]);
    fprintf('\n');
    
    disp(['The AUC from Z_AT: ' num2str(mean(auc_AD))]);
    disp(['The AUPR from Z_AT: ' num2str(mean(aupr_AD))]);
    
    fprintf('\n');
    disp(['The AUC from Z_A: ' num2str(mean(auc_A))]);
    disp(['The AUPR from Z_A: ' num2str(mean(aupr_A))]);
    
    fprintf('\n');
    disp(['The AUC from Z_ADT: ' num2str(mean(auc_ADT))]);
    disp(['The AUPR from Z_ADT: ' num2str(mean(aupr_ADT))]);
    
end

