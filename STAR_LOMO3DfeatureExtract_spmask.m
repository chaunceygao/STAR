function features = STAR_LOMO3DfeatureExtract_spmask(myopt, cam_num)
%% function STAR_LOMO3DfeatureExtract_spmask(myopt,cam_num)
% Function for LOMO3D feature extraction, by temperally aligned pooling
%
% Input:
%   <myopt>: parameters
%   <cam_num>: index of camera view
%
% Output:
%   <features>: LOMO3D features
%
% Date: 2016-07-18
% Author: Changxin Gao
% Institute: School of Automation, Huazhong University of Science and Technology
% Email: cgao@hust.edu.cn

path = myopt.path;
maskpath = myopt.maskpath;
cycle_flag = myopt.path;
len_cam = myopt.len_cam;
str = myopt.str;
camName = myopt.camName{cam_num}; 
path = [path camName];
maskpath = [maskpath camName];  

personpath = dir(path);

% load clip 
load(['clip\seq_video_clip_' camName(1:end-1) '_' str 'cycle.mat' ]);  % time_info

features = [];
for j=1:size(walking_cycle,1)
    person_index = cycle_index(j);
    imgpath = [path personpath(person_index+2).name];
    maskimgpath  = [maskpath personpath(person_index+2).name];
    imgfile = dir([imgpath '\*.png']);
    
    
    %% featrue extraction 
    for i=walking_cycle(j,1) : min(walking_cycle(j,end),length(imgfile))
        k = i - walking_cycle(j,1) + 1;
        images(:,:,:,k) = imread([imgpath '\' imgfile(i).name]);   
        
        msk = imread([maskimgpath '\' imgfile(i).name]);        
        masks(:,:,k) = msk;
    end  
    feature = STAR_LOMO3D_sp_mask(images,masks,myopt);
    features =  [features; feature'];   
end
