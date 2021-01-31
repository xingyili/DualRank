% Obtain_StatisticScore() -  compute the t statistics of genes
%
% Usage:
%              >> Obtain_StatisticScore(Genes_expression_value_cluster, Label_cluster, UniqueNode,Genes_ID);
%              >> [StatisticScore,stats] = ...
%                    Obtain_StatisticScore(Genes_expression_value_cluster, Label_cluster, UniqueNode,Genes_ID);
% Inputs:
%       Genes_expression_value_cluster   = Gene expression data of samples in each cluster
%       Label_cluster   = Labels of samples in each cluster
%       UniqueNode   = Genes in the network
%       Genes_ID   = Genes in the gene expression profile
% Outputs:
%       StatisticScore   = Statistic scores of genes
function StatisticScore = Obtain_StatisticScore(Genes_expression_value_cluster, Label_cluster, UniqueNode,Genes_ID)
    Genes_expression_data = Genes_expression_value_cluster;
    Label= Label_cluster;   
    [~,n1] = find(Label==0); %Good
    [~,n2] = find(Label==1); %Poor
    GoodSamples = Genes_expression_data(:,n1);
    PoorSamples = Genes_expression_data(:,n2); 
    StatisticScore = zeros(size(UniqueNode,1),2);
    StatisticScore(:,1) = UniqueNode;
    for k = 1:size(Genes_ID,1)
        [~,~,~,stats]=ttest2(GoodSamples(k,:),PoorSamples(k,:));  
        [~,~,kUniqueNode] = intersect(Genes_ID(k,1),UniqueNode);
        StatisticScore(kUniqueNode,2) = abs(stats.tstat);   
    end 
end