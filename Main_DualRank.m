clear
clc
dataPath = './data/';
[GenesExpressionData,~] = xlsread([dataPath,'GeneExpressionData.xlsx']);    
Genes_ID = GenesExpressionData(2:end,1);
label = GenesExpressionData(1,2:end);  
Genes_expression_value = GenesExpressionData(2:end,2:end);  

load([dataPath,'SparseAdjMatrix.mat'])
AdjSet.SparseAdjMatrix1 = SparseAdjMatrix1;     %AdjSet:adjacency matrix of networks
AdjSet.SparseAdjMatrix2 = SparseAdjMatrix2;
AdjSet.SparseAdjMatrix3 = SparseAdjMatrix3;
AdjSet.SparseAdjMatrix4 = SparseAdjMatrix4;
AdjSet.SparseAdjMatrix5 = SparseAdjMatrix5;
AdjSet.SparseAdjMatrix6 = SparseAdjMatrix6;

global Global_Var_RWRM
Global_Var_RWRM = [];

rs = 0.7;                                       % the parameter that quantifies the relative importance of the initial state
pro_jump = 0.5;                                 % the parameter that denotes the relative importance of the states of the neighbors in the layer or the states of nodes in other layers

%% ****************************PCA and GMM*****************************************
[clusterX,Gene_expression_data_normalized,gm_c1] = Obtain_ClusterSample(Genes_expression_value);  
%% ***************************Network propagation**********************************
NumberOfSample = size(clusterX,1);
Genes_expression_value_cluster1 = [];
Genes_expression_value_cluster2= [];
Genes_expression_value_cluster3= [];
Genes_expression_value_cluster4= [];
Label_cluster1= [];
Label_cluster2= [];
Label_cluster3= [];
Label_cluster4= [];
%%%%%%%%%%%%%%%%%%%%Gene expression and labels of samples in each cluster %%%%%%%%%%%%%%%%%
for i = 1:NumberOfSample
    if clusterX(i,1) == 1
       Genes_expression_value_cluster1_temp = Gene_expression_data_normalized(:,i);
       Label_cluster1_temp = label(1,i);
       Genes_expression_value_cluster1 = [Genes_expression_value_cluster1,Genes_expression_value_cluster1_temp];
       Label_cluster1 = [Label_cluster1,Label_cluster1_temp];
    elseif clusterX(i,1) == 2
       Genes_expression_value_cluster2_temp = Gene_expression_data_normalized(:,i);
       Label_cluster2_temp = label(1,i);
       Genes_expression_value_cluster2 = [Genes_expression_value_cluster2,Genes_expression_value_cluster2_temp];
       Label_cluster2 = [Label_cluster2,Label_cluster2_temp];
    elseif clusterX(i,1) == 3
       Genes_expression_value_cluster3_temp = Gene_expression_data_normalized(:,i);
       Label_cluster3_temp = label(1,i);
       Genes_expression_value_cluster3 = [Genes_expression_value_cluster3,Genes_expression_value_cluster3_temp];
       Label_cluster3 = [Label_cluster3,Label_cluster3_temp];
    else
       Genes_expression_value_cluster4_temp = Gene_expression_data_normalized(:,i);
       Label_cluster4_temp = label(1,i);
       Genes_expression_value_cluster4 = [Genes_expression_value_cluster4,Genes_expression_value_cluster4_temp];
       Label_cluster4 = [Label_cluster4,Label_cluster4_temp];       
    end
end
%%%%%%%%%%%%%%%%%%%%Network propagation%%%%%%%%%%%%%%%%%
for ClusterNumber = 1:gm_c1                                                %  gm_c1<=4
    if ClusterNumber == 1
        StatisticScore = Obtain_StatisticScore(Genes_expression_value_cluster1, Label_cluster1, UniqueNode,Genes_ID); 
        TableScores1 = A_RWR_M(AdjSet, StatisticScore(:,2), rs, pro_jump, false, [], 'row');
        Rank_Statistic = TableScores1.ScoreGeoMean;
        RankStatistic_1 = zeros(size(UniqueNode,1),2); 
        RankStatistic_1(:,1) = UniqueNode;
        RankStatistic_1(:,2) = Rank_Statistic;        
        Score_AllGene = RankStatistic_1(:,2) ;
        ScoreOfAllGene_1 = zeros(size(UniqueNode,1),2);  
        ScoreOfAllGene_1(:,1) = UniqueNode;   
        ScoreOfAllGene_1(:,2) = Score_AllGene; 
    elseif ClusterNumber == 2
        StatisticScore = Obtain_StatisticScore(Genes_expression_value_cluster2, Label_cluster2, UniqueNode,Genes_ID);
        TableScores2 = A_RWR_M(AdjSet, StatisticScore(:,2), rs,  pro_jump, true, [], 'row');
        Rank_Statistic = TableScores2.ScoreGeoMean;
        RankStatistic_2 = zeros(size(UniqueNode,1),2); 
        RankStatistic_2(:,1) = UniqueNode;
        RankStatistic_2(:,2) = Rank_Statistic;        
        Score_AllGene = RankStatistic_2(:,2) ;
        ScoreOfAllGene_2 = zeros(size(UniqueNode,1),2);  
        ScoreOfAllGene_2(:,1) = UniqueNode;   
        ScoreOfAllGene_2(:,2) = Score_AllGene;   
    elseif ClusterNumber == 3
        StatisticScore = Obtain_StatisticScore(Genes_expression_value_cluster3, Label_cluster3, UniqueNode,Genes_ID);
        TableScores3 = A_RWR_M(AdjSet, StatisticScore(:,2), rs,  pro_jump, true, [], 'row');
        Rank_Statistic = TableScores3.ScoreGeoMean;
        RankStatistic_3 = zeros(size(UniqueNode,1),2); 
        RankStatistic_3(:,1) = UniqueNode;
        RankStatistic_3(:,2) = Rank_Statistic;        
        Score_AllGene = RankStatistic_3(:,2);
        ScoreOfAllGene_3 = zeros(size(UniqueNode,1),2);  
        ScoreOfAllGene_3(:,1) = UniqueNode;   
        ScoreOfAllGene_3(:,2) = Score_AllGene;    
    else
        StatisticScore = Obtain_StatisticScore(Genes_expression_value_cluster4, Label_cluster4, UniqueNode,Genes_ID);    
        TableScores4 = A_RWR_M(AdjSet, StatisticScore(:,2), rs,  pro_jump, true, [], 'row');
        Rank_Statistic = TableScores4.ScoreGeoMean;
        RankStatistic_4 = zeros(size(UniqueNode,1),2); 
        RankStatistic_4(:,1) = UniqueNode;
        RankStatistic_4(:,2) = Rank_Statistic;        
        Score_AllGene = RankStatistic_4(:,2) ;
        ScoreOfAllGene_4 = zeros(size(UniqueNode,1),2);  
        ScoreOfAllGene_4(:,1) = UniqueNode;   
        ScoreOfAllGene_4(:,2) = Score_AllGene;          
    end
end
%% *******************Ranking*******************************************
NumOfFeature = 100;
clusternumber = gm_c1;
if clusternumber == 1
    [SortGene_1,~]=sortrows(ScoreOfAllGene_1,-2);
    Feature = unique(SortGene_1,'stable');
    [~,~,iGenes_ID] = intersect(Feature(1:NumOfFeature,:),Genes_ID);  
    GenesExpression = Genes_expression_value(iGenes_ID,:)';   
    CandidateFeature = Genes_ID(iGenes_ID,:)';
elseif clusternumber == 2
    [SortGene_1,~]=sortrows(ScoreOfAllGene_1,-2);
    [SortGene_2,~]=sortrows(ScoreOfAllGene_2,-2);
    Num = size(ScoreOfAllGene_1,1);
    Sort = zeros(Num,2);
    SortTemp = [];
    for kkk = 1:Num
        [~,~,iSortGene_2] = intersect(SortGene_1(kkk,1),SortGene_2(:,1));
        RankScore = (kkk+iSortGene_2)/2;
        SortTemp = [SortTemp;RankScore];
    end
    Sort(:,1) = SortGene_1(:,1);
    Sort(:,2) = SortTemp;
    SortNew = sortrows(Sort,2);
    Feature = SortNew(:,1);
    [~,~,iGenes_ID] = intersect(Feature(1:NumOfFeature,:),Genes_ID);  
    GenesExpression = Genes_expression_value(iGenes_ID,:)';
    CandidateFeature = Genes_ID(iGenes_ID,:)';
elseif clusternumber == 3
    [SortGene_1,~]=sortrows(ScoreOfAllGene_1,-2);
    [SortGene_2,~]=sortrows(ScoreOfAllGene_2,-2);
    [SortGene_3,~]=sortrows(ScoreOfAllGene_3,-2);
    Num = size(ScoreOfAllGene_1,1);
    Sort = zeros(Num,2);
    SortTemp = [];
    for kkk = 1:Num
        [~,~,iSortGene_2] = intersect(SortGene_1(kkk,1),SortGene_2(:,1));
        [~,~,iSortGene_3] = intersect(SortGene_1(kkk,1),SortGene_3(:,1));
        RankScore = (kkk+iSortGene_2+iSortGene_3)/3;
        SortTemp = [SortTemp;RankScore];
    end
    Sort(:,1) = SortGene_1(:,1);
    Sort(:,2) = SortTemp;
    SortNew = sortrows(Sort,2);
    Feature = SortNew(:,1);
    [~,~,iGenes_ID] = intersect(Feature(1:NumOfFeature,:),Genes_ID); 
    GenesExpression = Genes_expression_value(iGenes_ID,:)';
    CandidateFeature = Genes_ID(iGenes_ID,:)';
else
    [SortGene_1,~]=sortrows(ScoreOfAllGene_1,-2);
    [SortGene_2,~]=sortrows(ScoreOfAllGene_2,-2);
    [SortGene_3,~]=sortrows(ScoreOfAllGene_3,-2);
    [SortGene_4,~]=sortrows(ScoreOfAllGene_4,-2);
    Num = size(ScoreOfAllGene_1,1);
    Sort = zeros(Num,2);
    SortTemp = [];
    for kkk = 1:Num
        [~,~,iSortGene_2] = intersect(SortGene_1(kkk,1),SortGene_2(:,1));
        [~,~,iSortGene_3] = intersect(SortGene_1(kkk,1),SortGene_3(:,1));
        [~,~,iSortGene_4] = intersect(SortGene_1(kkk,1),SortGene_4(:,1));
        RankScore = (kkk+iSortGene_2+iSortGene_3+iSortGene_4)/4;
        SortTemp = [SortTemp;RankScore];
    end
    Sort(:,1) = SortGene_1(:,1);
    Sort(:,2) = SortTemp;
    SortNew = sortrows(Sort,2);
    Feature = SortNew(:,1);
    [~,~,iGenes_ID] = intersect(Feature(1:NumOfFeature,:),Genes_ID); 
    GenesExpression = Genes_expression_value(iGenes_ID,:)';
    CandidateFeature = Genes_ID(iGenes_ID,:)';
end
LLabel = label'; 
group = (LLabel==0); 
%% *****************************AUC optimization model****************************
feature_index = NCC_AUC(GenesExpression,group,0.00001);   % indexes of biomarkers     
Biomarker = CandidateFeature(:,feature_index)'; 
save Biomarker.mat Biomarker
