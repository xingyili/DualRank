% Obtain_ClusterSample() - PCA and sample cluster using GMM 
%
% Usage:
%              >> Obtain_ClusterSample(Training_ExpressionData);
%              >> [clusterX,Gene_expression_data_normalized,gm_c1] = ...
%                    Obtain_ClusterSample(Training_ExpressionData);
% Inputs:
%       Training_ExpressionData   = Gene expression data of training samples
% Outputs:
%       clusterX   = cluster of each sample
%       Gene_expression_data_normalized = Gene expression data after normalization
%       gm_c1  = Number of cluster
function [clusterX,Gene_expression_data_normalized,gm_c1] = Obtain_ClusterSample(Training_ExpressionData)
    [Gene_expression_data_normalized,~,~] = zscore(Training_ExpressionData);
    Gene_expression_data_normalized_trans = Gene_expression_data_normalized';
    [~,SCORE,~] = pca(Gene_expression_data_normalized_trans);  
    pcaData1 = SCORE(:,1:2);           

    k = 1:3;
    nK = numel(k);
    Sigma = {'diagonal','full'};
    nSigma = numel(Sigma);
    SharedCovariance = {true,false};
    SCtext = {'true','false'};
    nSC = numel(SharedCovariance);
    RegularizationValue = 0.01;
    options = statset('MaxIter',10000);

    gm = cell(nK,nSigma,nSC);
    aic = zeros(nK,nSigma,nSC);
    bic = zeros(nK,nSigma,nSC);
    converged = false(nK,nSigma,nSC);

    for m = 1:nSC
        for j = 1:nSigma
            for i = 1:nK
                gm{i,j,m} = fitgmdist(pcaData1,k(i),...
                    'CovarianceType',Sigma{j},...
                    'SharedCovariance',SharedCovariance{m},...
                    'RegularizationValue',RegularizationValue,...
                    'Options',options);
                aic(i,j,m) = gm{i,j,m}.AIC;
                bic(i,j,m) = gm{i,j,m}.BIC;
                converged(i,j,m) = gm{i,j,m}.Converged;
            end
        end
    end

    aic_2D = reshape(aic,nK,nSigma*nSC);
    bic_2D = reshape(bic,nK,nSigma*nSC);
    aic_bic_sum = aic_2D + bic_2D;
    [row_a,col_b] = find(aic_bic_sum==min(min(aic_bic_sum)));
    gm_c1 = row_a(1,1);
    col_b = col_b(1,1);
    if col_b<=size(aic_bic_sum,2)/2
        gm_c2 = col_b;
        gm_c3 = 1;
    else
        gm_c2 = col_b - size(aic_bic_sum,2)/2;
        gm_c3 = 2;
    end
    gmBest = gm{gm_c1,gm_c2,gm_c3};
    clusterX = cluster(gmBest,pcaData1);
end