function [test, train_sum_hist, updata_num, update_sp_num, update_flag, update_negtive]... 
          = t1_update_app_model(myopt, update, test, update_hist_sum, update_index_pre, update_index_pre_final)
%% Copyright (C) Shu Wang.
%% All rights reserved.

test.update_interval_num = 0;
test.save_prob = sum(test.update_spt_conf) / myopt.update_incre_num;
test.save_std = std(test.update_spt_conf(:),1);

clear test.clust_Cent test.data2cluster test.cluster2dataCell test.TrainCluster_Weight;

%% superpixels clustering using Mean-Shift
[test.clust_Cent, test.data2cluster, test.cluster2dataCell] = MeanShiftCluster(update_hist_sum, myopt.cluster_bandWidth);        
cluster_Sum = size(test.cluster2dataCell, 1);      % cluster number
test.train_sp_sum_num = update_index_pre(myopt.update_incre_num);

test.TrainCluster_Weight = t1_show_cal_cluster_wt_update(cluster_Sum, test.data2cluster, myopt, update, update_index_pre_final);
train_sum_hist = update_hist_sum;
updata_num = 0;
update_sp_num = 0;
update_flag = 1;
update_negtive = 0;

