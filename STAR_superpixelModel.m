function [person] = STAR_superpixelModel(myopt,cam_num)
%% function STAR_superpixelModel(myopt,cam_num)
% Function for Motion information extraction, by selecting superpixels for tracking
%
% Input:
%   <myopt>: parameters
%   <cam_num>: index of camera view
%
% Output:
%   <person>: candidate motion trajectory of superpixels
%
% Date: 2016-07-18
% Author: Changxin Gao
% Institute: School of Automation, Huazhong University of Science and Technology
% Email: cgao@hust.edu.cn

%% set paths for images and masks
path = myopt.path;
maskpath = myopt.maskpath;
camName = myopt.camName{cam_num};   
path = [path camName];
maskpath = [maskpath  camName]; 

personpath = dir(path);

% find the videos, whose frame number is bigger than a given threshold (the frame nubmers of a cycle). 
% If the frame number of a video is smaller than the threshold, all the frames are used.
len_cam = myopt.len_cam;
personlist = len_cam(:,1)>myopt.frameNumber & len_cam(:,2)>myopt.frameNumber;

for j=1:size(personlist,1)%length(personpath)-2
    imgpath = [path personpath(j+2).name];
    maskimgpath = [maskpath personpath(j+2).name];
    imgfile = dir([imgpath '\*.png']);
    
    % if the frame number of a video is smaller than the threshold 
    if personlist(j) == 0
        person(j).posDiff = [];
        continue;
    end   
    
    train_sp_sum_num = 0;    
    sp_index_pre(1) = 0;
    
    seq = [];
    %% (1) superpixel segmentation   and  color  featrue extraction   
    for i=1:length(imgfile)
        seq(i).warpimg_rgb = imread([imgpath '\' imgfile(i).name]);
        seq(i).mask = imread([maskimgpath '\' imgfile(i).name]);
        
        % superpixel segmentation
        seq(i).labels = SLIC_mex(seq(i).warpimg_rgb, myopt.SLIC_sp_num, myopt.SLIC_spatial_proximity_weight);
        
        N_superpixels = unique(seq(i).labels);
        seq(i).sp_num = max(N_superpixels);               
        
        % rows and cols of the superpixels
        for k=1:seq(i).sp_num
            [row col] = find(seq(i).labels==k);
            seq(i).row(k) = mean(row);
            seq(i).col(k) = mean(col);
        end
        
        % calculate superPixel histogram
        train_sp_sum_num = train_sp_sum_num + seq(i).sp_num;
        sp_index_pre(i+1) = train_sp_sum_num;
        
        % calculate the HSV color histogram of this frame 
        seq(i).warpimg_hsi = rgb2hsi(seq(i).warpimg_rgb);
        [~, temp_index, temp_sp_cl_hist] = t1_cal_hsi_hist(seq(i), i, myopt.ch_bins_num, N_superpixels);
        seq(i).sp_cl_hist = temp_sp_cl_hist;      % the histogram of this frame        
    end       

     
     %% tracking selected superpixels to extract motion information     
      posDiff = [];
      for i=1:1%length(seq) - 30  % start from 1 to end-60
%           if size(posDiff,1) > 30
%               continue;
%           end
          % A superpixel is considered to track, if its overlap with mask
          % is larger than 0.5, and it is in the lower part
          for p = 1:seq(i).sp_num
              
              indNum = p;
%               temp_pixel_ratio_mask = sum(seq(i).mask(find(seq(i).labels==indNum))>0) / length(find(seq(i).labels==indNum));
              temp_pixel_ratio_mask = 1;

              if seq(i).row(indNum)<100  |  temp_pixel_ratio_mask < 0.5
                  continue;
              end              
              sp_seq_pos_row = zeros(1,length(seq)) + myopt.row/2;
              sp_seq_pos_col = zeros(1,length(seq)) + myopt.col/2;
              sp_seq_pos_row(i) = seq(i).row(indNum);
              sp_seq_pos_col(i) = seq(i).col(indNum);

              for k = i+1:length(seq)
                  % Method(1)  template matching
                  tracking_ind = i;
                  temp_index = indNum;                  
                  % Mehtod(2)  tracking
                  % tracking_ind = k-1;
                  
                  % distances to the superpixels
                  dis_to_others = slmetric_pw(seq(tracking_ind).sp_cl_hist(:,temp_index), seq(k).sp_cl_hist, 'sqdist');
                  pos_dis_row = abs(sp_seq_pos_row(k-1) - seq(k).row); % posistion distance
                  pos_dis_col = abs(sp_seq_pos_col(k-1) - seq(k).col);         
                  
                  mask_dis = zeros(1,length(pos_dis_col));
                  for q=1:length(pos_dis_row)
                      mask_dis(q) = seq(k).mask(uint16(seq(k).row(q)), uint16(seq(k).col(q)));
                  end

                  disp([j i p k])                  
                  % the factors considered to find the best match
                  combine_dis =  dis_to_others + double(pos_dis_row>5) + double(1-mask_dis) ...
                      + double(pos_dis_col>15) ;% .* exp(-pos_dis_col);
                  
                  [~, temp_index] = min(combine_dis, [], 2);
                  sp_seq_pos_row(k) = seq(k).row(temp_index);
                  sp_seq_pos_col(k) = seq(k).col(temp_index);                 
 
              end
              posDiff = [posDiff; sp_seq_pos_col];              
          end        
      end
      person(j).posDiff = posDiff;
end