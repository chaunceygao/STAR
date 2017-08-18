function STAR()
%% This is a demo for the STAR method for re-id on the iLIDS-VID database. 
% This code is for the our submitted manuscript: 
% Changxin Gao, Jin Wang, Leyuan Liu, Jin-Gang Yu, Nong Sang. Superpixel-based 
% Temporally Aligned Representation for Video-based Person Re-identification
% 
% You can run this script to reproduce our results. 
% 
% Date: 2016-07-18
% Author: Changxin Gao
% Institute: School of Automation, Huazhong University of Science and Technology
% Email: cgao@hust.edu.cn

clear
close all

addpath(genpath('.\'))
run('.\utils\vlfeat\vlfeat-0.9.20\toolbox\vl_setup.m')

%% parameters initial
myopt = STAR_parameterInitial();
str = myopt.str;

%% Person segmentation
for cam_num = 1:2  % camera num
    STAR_personSegmentation( myopt, cam_num );
end

%% Motion information extraction, by selecting superpixels for tracking
for cam_num = 1:2  % camera num
    person = STAR_superpixelModel(myopt,cam_num);
    save(['positions\pos_' myopt.camName{cam_num}(1:end-1) '.mat' ],'person');
end


%% walking cycle selection, e.g., time alienment by motion information (position (x) curves) processing
for cam_num = 1:2  % camera num
    [walking_cycle cycle_index]= STAR_timealienment(myopt, cam_num);%
    save(['clip\seq_video_clip_' myopt.camName{cam_num}(1:end-1) '_' str 'cycle.mat' ],'walking_cycle', 'cycle_index');
end

%% Temporally Aligned Superpixel-based Representation (STAR) for videos
for cam_num = 1:2  % camera num
    features = STAR_LOMO3DfeatureExtract_spmask(myopt,cam_num);
    save(['features\features_' myopt.camName{cam_num}(1:end-1) '_' str 'cycle_STAR_' myopt.poolingmethod int2str(myopt.partnum) '.mat' ],'features', 'cycle_index');
end


%% performance evaluation
[meanCms] = STAR_performance(myopt) ;
save(['meanCms_STAR_' myopt.poolingmethod int2str(myopt.partnum) '.mat' ],'meanCms');