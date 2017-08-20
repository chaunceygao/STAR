function test = t1_construct_appearance_model(train_sum_hist, myopt, frame, train_sp_sum_num, sp_index_pre, f)
%% Copyright (C) Shu Wang.
%% All rights reserved.

%% superpixels clustering using Mean-Shift
[test.clust_Cent, test.data2cluster, test.cluster2dataCell] = MeanShiftCluster(train_sum_hist, myopt.cluster_bandWidth);      
test.cluster_Sum = size(test.cluster2dataCell, 1);      % number of clusters

%% calculate the confidence for each cluster to construct appearance model
test.TrainCluster_Weight = t1_show_cal_cluster_wt(test.cluster_Sum, test.data2cluster, myopt,frame, sp_index_pre);

%% parameters for the next tracking frame
test.block_size = [frame(f).p(4),frame(f).p(3)];
test.p = frame(f).p;    
test.est = affparam2ultimate(test.p ,test.block_size);
test.warpimg =  frame(f).warpimg;
test.warp_p = frame(f).warp_p; 
test.train_sp_sum_num = train_sp_sum_num;

