function [test, update, updata_num, update_index_pre, update_hist_sum, update_sp_num, update_index_pre_final] ...
          = t1_update_info(myopt, test, update, updata_num, update_index_pre, update_hist_sum, update_sp_num, update_index_pre_final)   
%% Copyright (C) Shu Wang.
%% All rights reserved.

%% accumulate information to update the appearance model
if test.update_flag == 1
    test.update_interval_num = test.update_interval_num+1;
    if mod(test.update_interval_num, myopt.update_spacing) == 0
        updata_num = updata_num + 1;
        test.incre_prob = test.incre_prob + test.spt_conf;
        if length(update_index_pre) < myopt.update_incre_num
            update(updata_num).labels = test.labels;
            update_hist_sum = [update_hist_sum, test.temp_sp_cl_hist];
            update(updata_num).warp_p = test.warp_p;
            update(updata_num).warpimg_tmpl = test.warpimg_tmpl;
            update_sp_num = update_sp_num + test.sp_num;
            update_index_pre = [update_index_pre, update_sp_num];
            update_index_pre_final = [0, update_index_pre];
            test.update_spt_conf = [test.update_spt_conf, test.spt_conf];
        else
            update(1) = [];
            temp_update = update;
            clear update;
            update = temp_update;
            clear temp_update;
            update_hist_sum = update_hist_sum(:,update_index_pre(1)+1:update_index_pre(myopt.update_incre_num));
            clear temp_hist;
            temp_hist = update_hist_sum;
            clear update_hist_sum;
            update_hist_sum = temp_hist;
            clear temp_hist;
            update_hist_sum = [update_hist_sum , test.temp_sp_cl_hist];

            update_index_pre = [update_index_pre, test.sp_num + update_index_pre(myopt.update_incre_num)];
            update_index_pre = update_index_pre - update_index_pre(1);
            update_index_pre_final = update_index_pre;
            update_index_pre(1)=[];
            temp_update = update_index_pre;
            clear update_index_pre;
            update_index_pre = temp_update;
            clear temp_update;                

            for i = 1:myopt.update_incre_num-1
                test.update_spt_conf(i) = test.update_spt_conf(i+1);
            end
            test.update_spt_conf(myopt.update_incre_num) = test.spt_conf;

            update(myopt.update_incre_num).labels = test.labels;
            update(myopt.update_incre_num).warp_p = test.warp_p;
            update(myopt.update_incre_num).warpimg_tmpl = test.warpimg_tmpl;
            update(myopt.update_incre_num).spt_conf = test.spt_conf;
        end

    end
else
    test.update_interval_num = test.update_interval_num+1;
    text(5, 50, ['Severe Occlusion!!'], 'Color','y', 'FontWeight','bold', 'FontSize',24);
    if mod(test.update_interval_num, myopt.update_spacing) == 0
        updata_num = updata_num + 1;
        test.incre_prob = test.incre_prob + test.save_prob;
        if length(update_index_pre) < myopt.update_incre_num
            update(updata_num).labels = test.labels;
            update_hist_sum = [update_hist_sum , test.temp_sp_cl_hist];
            update(updata_num).warp_p = test.warp_p;
            update(updata_num).warpimg_tmpl = test.warpimg_tmpl;
            update_sp_num = update_sp_num + test.sp_num;
            update_index_pre = [ update_index_pre, update_sp_num];
            update_index_pre_final = [0 , update_index_pre];
            test.update_spt_conf = [test.update_spt_conf, test.save_prob];
        else
            update(myopt.update_incre_num - 2) = [];
            update_hist_sum(:, update_index_pre(myopt.update_incre_num - 3)+1 : update_index_pre(myopt.update_incre_num -2)) = [];
            update_hist_sum = [update_hist_sum , test.temp_sp_cl_hist];

            update_index_pre(myopt.update_incre_num - 2 : myopt.update_incre_num ) = ...
            update_index_pre(myopt.update_incre_num - 2 : myopt.update_incre_num ) - (update_index_pre(myopt.update_incre_num - 2) ...
                                                                                   - update_index_pre(myopt.update_incre_num - 3));
            update_index_pre(myopt.update_incre_num - 2) = [];
            update_index_pre = [update_index_pre, test.sp_num + update_index_pre( myopt.update_incre_num-1)];
            update_index_pre_final = [0, update_index_pre];

            test.update_spt_conf(1) = [];
            test.update_spt_conf = [test.update_spt_conf, test.save_prob];

            update(myopt.update_incre_num).labels = test.labels;
            update(myopt.update_incre_num).warp_p = test.warp_p;
            update(myopt.update_incre_num).warpimg_tmpl = test.warpimg_tmpl;
            update(myopt.update_incre_num).spt_conf = test.spt_conf;
        end
        update(updata_num).warp_p(3) = 0;
        update(updata_num).warp_p(4) = 0;
   end
end 

