function [frame, train_sp_sum_num, sp_index_pre, train_sum_hist_index, train_sum_hist] ...
          = t1_train_info(opt, myopt, f, train_box_param, frame, train_sp_sum_num, ...
                          sp_index_pre, train_sum_hist_index, train_sum_hist)
%% Copyright (C) Shu Wang.
%% All rights reserved.

frame(f).p = affparam2original(train_box_param(:,f),opt.tmplsize); % 1*5 matrix 
temp_length = uint16(norm([frame(f).p(3)/2,frame(f).p(4)/2]) * myopt.grid_ratio);

frame(f).train_box = zeros(1,4);
frame(f).train_box(1) = max(1,frame(f).p(1) - temp_length);
frame(f).train_box(2) = max(1,frame(f).p(2) - temp_length);
frame(f).train_box(3) = min(myopt.image_size.cx , frame(f).p(1) + temp_length);
frame(f).train_box(4) = min(myopt.image_size.cy , frame(f).p(2) + temp_length);  

frame(f).warpimg =  frame(f).image(frame(f).train_box(2) : frame(f).train_box(4),...
                               frame(f).train_box(1) : frame(f).train_box(3), :);
warpimg_size = size( frame(f).warpimg);
frame(f).warpimg_tmpl.cx = warpimg_size(2);
frame(f).warpimg_tmpl.cy = warpimg_size(1);
frame(f).warp_p = frame(f).p;
frame(f).warp_p(1) = frame(f).warp_p(1) - frame(f).train_box(1);
frame(f).warp_p(2) = frame(f).warp_p(2) - frame(f).train_box(2);
frame(f).warp_param.est = affparam2ultimate(frame(f).warp_p , opt.tmplsize); 

%% segment the frames and give the labels
% SLIC segmentation  
frame(f).labels = SLIC_mex(frame(f).warpimg, myopt.SLIC_sp_num, myopt.SLIC_spatial_proximity_weight);

% frame(f).labels = readDAT(warpimg_size,'temp_image.dat');
% temp_I_s = imread( 'temp_image_SLIC.bmp');
% frame(f).I_s =  temp_I_s;
% for i = 1:warpimg_size(1)
%     for  j = 1:warpimg_size(2)
%         frame(f).I_s(i, j,:) = temp_I_s(i, (warpimg_size(2)+1 - j), :);
%     end
% end

N_superpixels = unique(frame(f).labels);
frame(f).sp_num = max(N_superpixels);       % record the number of superpixels of this frame

%% calculate superPixel histogram
train_sp_sum_num = train_sp_sum_num + frame(f).sp_num;
sp_index_pre(f+1) = train_sp_sum_num;

% calculate the HSI color histogram of this frame 
frame(f).warpimg_hsi = rgb2hsi(frame(f).warpimg);
[~, temp_index, temp_sp_cl_hist] = t1_cal_hsi_hist(frame(f), f, myopt.ch_bins_num, N_superpixels);
train_sum_hist_index = [train_sum_hist_index, temp_index]; % histo2pixel index        
frame(f).sp_cl_hist = temp_sp_cl_hist;      % the histogram of this frame
train_sum_hist = [train_sum_hist, frame(f).sp_cl_hist];   % the histogram of all training frames    

