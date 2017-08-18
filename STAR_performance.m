function [meanCms] = STAR_performance(myopt) 
%% function STAR_performance(myopt)
% Function for performance evaluation 
% We evaluate the performance using the code of the paper: 
%   Shengcai Liao, Yang Hu, Xiangyu Zhu, and Stan Z. Li. Person
%   re-identification by local maximal occurrence representation and metric
%   learning. In IEEE Conference on Computer Vision and Pattern Recognition, 2015.
%
% Date: 2016-07-18
% Author: Changxin Gao
% Institute: School of Automation, Huazhong University of Science and Technology
% Email: cgao@hust.edu.cn

numClass = myopt.numClass;
numFolds = myopt.numFolds;
numRanks = myopt.numRanks;

cycle_flag = myopt.cycle_flag;
ls_set = myopt.ls_set;

if cycle_flag == 1
    str = 'one';
elseif cycle_flag == 0.5
    str = 'half';
end

%% load the features
outF = importdata(['features\features_' myopt.camName{1}(1:end-1) '_' str 'cycle_STAR_' myopt.poolingmethod int2str(myopt.partnum) '.mat'] );
galFea = outF.features;
galLabel = outF.cycle_index;
outF = importdata(['features\features_' myopt.camName{2}(1:end-1) '_' str 'cycle_STAR_' myopt.poolingmethod int2str(myopt.partnum) '.mat']);
probFea = outF.features;
probLabel = outF.cycle_index;
clear outF

%% evaluate
cms = zeros(numFolds, numRanks);

for nf = 1 : numFolds
    p = ls_set(nf,:);
    
    galFea1 = galFea( p(1:numClass/2), : );
    probFea1 = probFea( p(1:numClass/2), : );
    
    t0 = tic;
    [W, M] = XQDA(galFea1, probFea1, (1:numClass/2)', (1:numClass/2)');
%     options.qdaDims = 100;
%     [W, M] = XQDA(galFea1, probFea1, (1:numClass/2)', (1:numClass/2)', options);
    
    %{
    %% if you need to set different parameters other than the defaults, set them accordingly
    options.lambda = 0.001;
    options.qdaDims = -1;
    options.verbose = true;
    [W, M] = XQDA(galFea1, probFea1, (1:numClass/2)', (1:numClass/2)', options);
    %}
    
    clear galFea1 probFea1
    trainTime = toc(t0);
    
    galFea2 = galFea(p(numClass/2+1 : end), : );
    probFea2 = probFea(p(numClass/2+1 : end), : );
    
    t0 = tic;
    dist = MahDist(M, galFea2 * W, probFea2 * W);
    clear galFea2 probFea2 M W
    matchTime = toc(t0);
    
    fprintf('Fold %d: ', nf);
    fprintf('Training time: %.3g seconds. ', trainTime);    
    fprintf('Matching time: %.3g seconds.\n', matchTime);
    
    cms(nf,:) = EvalCMC( -dist, 1 : numClass / 2, 1 : numClass / 2, numRanks );
    clear dist
    
    fprintf(' Rank1,  Rank5, Rank10, Rank15, Rank20\n');
    fprintf('%5.2f%%, %5.2f%%, %5.2f%%, %5.2f%%, %5.2f%%\n\n', cms(nf,[1,5,10,15,20]) * 100);
end

meanCms = mean(cms);
plot(1 : numRanks, meanCms);

fprintf('The average performance:\n');
fprintf(' Rank1,  Rank5, Rank10, Rank15, Rank20\n');
fprintf('%5.2f%%, %5.2f%%, %5.2f%%, %5.2f%%, %5.2f%%\n\n', meanCms([1,5,10,15,20]) * 100);