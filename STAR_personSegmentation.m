function STAR_personSegmentation(myopt,cam_num)
%% function STAR_personSegmentation(myopt,cam_num)
% Function for the person segmentation (parsing) 
% We segment persons using the code (http://mmlab.ie.cuhk.edu.hk/projects/luoWTiccv2013DDN/index.html)
% of the paper: P. Luo, X. Wang, and X. Tang, Pedestrian Parsing via Deep Decompositional Neural Network, 
% in Proceedings of IEEE International Conference on Computer Vision (ICCV) 2013
%
% Input:
%   <myopt>: parameters
%   <myopt>: index of camera view
% 
% Date: 2016-07-18
% Author: Changxin Gao
% Institute: School of Automation, Huazhong University of Science and Technology
% Email: cgao@hust.edu.cn

%% load models
path = '.\PedParsing\'; % your path '.\PedParsing\'
load([ path 'model\coarse_compressed.mat']);
model_coarse.W = Uk*full(Sk)*Vk'; model_coarse.b = b;
load([ path 'model\fine_compressed.mat']);
model_fine.W = Uk*full(Sk)*Vk'; model_fine.b = b;
load([ path  'model\clrmap.mat']);

%% set paths for images and masks
path = myopt.path;
maskpath = myopt.maskpath;
camName = myopt.camName{cam_num};   
t_path = [path camName];
t_path_mask = [maskpath  camName]; 
mkdir(t_path_mask);
personpath = dir(t_path);

%% person segmentation (parsing)
for j=1:length(personpath)-2
    imgpath = [t_path personpath(j+2).name];
    imgpath_mask = [t_path_mask personpath(j+2).name];
    mkdir(imgpath_mask);
    
    imgfile = dir([imgpath '\*.png']);
    
    for i=1:length(imgfile)
        img = imread([imgpath '\' imgfile(i).name]);
        [h w d] = size(img);
        img = imresize(img,[160 60]);
        labelmap = decompositionalLayer(img, model_coarse, model_fine);
        imwrite(imresize(labelmap, 2, 'nearest'), colormap, [imgpath_mask '\' imgfile(i).name]);
        img = imread([imgpath_mask '\' imgfile(i).name]);
        img = imresize(img,[h w],'nearest');
        imwrite(img, [imgpath_mask '\' imgfile(i).name]);  
    end    
end