function TrainCluster_Weight = t1_show_cal_cluster_wt(cluster_Sum, data2cluster, myopt, frame, sp_index_pre)
%% Copyright (C) Shu Wang.
%% All rights reserved.
myopt.train_frame_num = min(myopt.num_frame_inf,length(frame));

train_frame_num = myopt.train_frame_num;
negative_penalty_ratio = myopt.negative_penalty_ratio;
TrainCluster_Weight = zeros(1,cluster_Sum);  
temp_train_cluster_wt = zeros(2,cluster_Sum);

for f = 1:train_frame_num
    temp_labels = frame(f).labels;
    tmpl = frame(f).warpimg_tmpl;
    temp_warp_p = frame(f).warp_p;
    costheta = cos(temp_warp_p(5));
    sintheta = sin(temp_warp_p(5));
    x0 = costheta * temp_warp_p(1) + sintheta * temp_warp_p(2);
    y0 = costheta * temp_warp_p(2) + sintheta * temp_warp_p(1);
    for i = 1:tmpl.cy     
        for j = 1:tmpl.cx  
            pixel_index = sp_index_pre(f) + temp_labels(i,j);
            if pixel_index  > length(data2cluster)
                pixel_index = length(data2cluster);
            end
            temp_cluster_num = data2cluster(pixel_index);
            if abs(costheta*j+sintheta*i - x0) < 0.5 * temp_warp_p(3) && abs(costheta*i+sintheta*j - y0) < 0.5 * temp_warp_p(4)
                % sum of positive pixels
                temp_train_cluster_wt(1,temp_cluster_num) = temp_train_cluster_wt(1,temp_cluster_num) + 1;
            else
                % sum of negative pixels
                temp_train_cluster_wt(2,temp_cluster_num) = temp_train_cluster_wt(2,temp_cluster_num) + 1;
            end
        end
    end
end

temp_train_cluster_wt = double(temp_train_cluster_wt);

for i = 1:cluster_Sum
    TrainCluster_Weight(i,1) = double((temp_train_cluster_wt(1,i) - temp_train_cluster_wt(2,i)) / (temp_train_cluster_wt(1,i) + temp_train_cluster_wt(2,i)));
    TrainCluster_Weight(i,2) = double( temp_train_cluster_wt(1,i) / (temp_train_cluster_wt(1,i) + temp_train_cluster_wt(2,i)));
    if TrainCluster_Weight(i,1) < 0
        TrainCluster_Weight(i,1) = TrainCluster_Weight(i,1)/negative_penalty_ratio;
    end
end
