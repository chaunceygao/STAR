function myopt = STAR_parameterInitial()
%% function STAR_parameterInitial() 
% Date: 2016-07-18
% Author: Changxin Gao
% Institute: School of Automation, Huazhong University of Science and Technology
% Email: cgao@hust.edu.cn

%% modify these parameters for different datasets
myopt.path = '.\iLIDS-VID\sequences\'; % change to YOUR PATH
myopt.maskpath = '.\iLIDS-VID\masks\'; % change to YOUR PATH
myopt.camName{1} = ['cam1\'];
myopt.camName{2} = ['cam2\'];
myopt.row = 128;
myopt.col = 64;

%% 1 for one cycle, 0.5 for half cycle
myopt.cycle_flag = 1;
if myopt.cycle_flag == 1
    myopt.str = 'one';
elseif myopt.cycle_flag == 0.5
    myopt.str = 'half';
end

%% parameters for person segmentation
myopt.stripe = 5 ;
myopt.regionSize = [3 3 3 3 3 3] ;
myopt.regularizer = 0.1 ;
myopt.SLIC_sp_num = 100 ;
myopt.SLIC_spatial_proximity_weight = 40 ;
myopt.ch_bins_num = 8;

%% parameters for time alignment
myopt.frameNumber = 26;
myopt.threshold = 15;

%% parameters for feature extraction
myopt.numScales = 3;
myopt.blockSize = 10;
myopt.blockStep = 5;

myopt.regionSize = [3 3 3 3 3 3] ;
myopt.regularizer = 0.1 ;

myopt.hsvBins = [8,8,8];
myopt.tau = 0.3;
myopt.R = [3, 5];
myopt.numPoints = 4;

%% parameters for temporally aligned pooling
% myopt.poolingmethod = 'max';
% myopt.poolingmethod = 'key';
myopt.poolingmethod = 'avg';
myopt.partnum = 8;

myopt.threshold = 15;

%% multiple walking cycles or one single cycle for each person
myopt.numcycle = 'single';
% myopt.numcycle = 'multiple';

%% training & testing
myopt.numClass = 300;
myopt.numFolds = 10;
myopt.numRanks = 100;
load('.\splits\train_test_splits_ilidsvid.mat');
myopt.ls_set = ls_set;

%% length of image sequences
myopt.showresult = 0;
% load('len_cam');
len_cam = [];
for cam_num = 1:2  % camera num
    camName = myopt.camName{cam_num}; 
    curpath = [myopt.path camName];
    personpath = dir(curpath);
    for j=1:length(personpath)-2
        imgpath = [curpath personpath(j+2).name];
        imgfile = dir([imgpath '\*.png']);
        len_cam(j,cam_num) = length(imgfile);
    end
end
myopt.len_cam = len_cam;