# DualRank

A dual ranking algorithm based on the multiplex network (DualRank) toward heterogeneous complex disease analysis

Main_DualRank.m is the main script to call DualRank.

There are some Matlab scripts for each step of DualRank analysis, and called in Main_DualRank.m:

  i. Obtain_ClusterSample.m 
  ii. Obtain_StatisticScore.m 
  iii. A_RWR_M.m 
  iv. A_RWRplus.m 
  v. NCC_AUC.m 
  vi. getNormalizedMatrix.m 
  vii. getMultiplexMatrixFromAdjSet.m 

The input datasets include:

  i. GeneExpressionData.xlsx 
  ii. SparseAdjMatrix.mat

The final biomarker list is saved as Biomarker.mat

As a demo, users can directly run Main_DualRank.m in Matlab. This package has been tested in different computer environments as: Window 7 or above; Matlab 2016 or above.
