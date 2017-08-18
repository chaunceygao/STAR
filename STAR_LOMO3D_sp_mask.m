function descriptors = STAR_LOMO3D_sp_mask(images, masks, myopt, options)
%% function Descriptors = STAR_LOMO3D_sp_mask(images, masks, myopt, options)
% Function for the Superpixel-based Temprally Aligned Representation (STAR) extraction
% We modify the code of the paper: 
%   Shengcai Liao, Yang Hu, Xiangyu Zhu, and Stan Z. Li. Person
%   re-identification by local maximal occurrence representation and metric
%   learning. In IEEE Conference on Computer Vision and Pattern Recognition, 2015.
%
% Input:
%   <images>: a set of n RGB color images. Size: [h, w, 3, n]
%   <masks>: a set of n masks. Size: [h, w, n]
%   <myopt>: parameters
%   [optioins]: optional parameters. A structure containing any of the
%   following fields:
%       numScales: number of pyramid scales in feature extraction. Default: 3
%       blockSize: size of the sub-window for histogram counting. Default: 10
%       blockStep: sliding step for the sub-windows. Default: 5
%       hsvBins: number of bins for HSV channels. Default: [8,8,8]
%       tau: the tau parameter in SILTP. Default: 0.3
%       R: the radius paramter in SILTP. Specify multiple values for multiscale SILTP. Default: [3, 5]
%       numPoints: number of neiborhood points for SILTP encoding. Default: 4
%   The above default parameters are good for 128x48 and 160x60 person
%   images. You may need to adjust the numScales, blockSize, and R parameters
%   for other smaller or higher resolutions.
%
% Output:
%   descriptors: the extracted STAR descriptors. Size: [d, n]
% 
% Date: 2016-07-18
% Author: Changxin Gao
% Institute: School of Automation, Huazhong University of Science and Technology
% Email: cgao@hust.edu.cn


poolingmethod = myopt.poolingmethod ;
partnum = myopt.partnum ;

numScales = myopt.numScales ;
blockSize = myopt.blockSize ;
blockStep = myopt.blockStep ;

hsvBins = myopt.hsvBins ;
tau = myopt.tau ;
R = myopt.R ;
numPoints = myopt.numPoints ;

regionSize = myopt.regionSize ;
regularizer = myopt.regularizer ;



timeSize = 1/partnum;  % 
timeStep = 1/partnum;  % 
timestart = 0:timeStep:1-timeSize; 
timeend = timestart + timeSize;
keynum = 0:timeStep:1;

if nargin >= 4
    if isfield(options,'numScales') && ~isempty(options.numScales) && isscalar(options.numScales) && isnumeric(options.numScales) && options.numScales > 0
        numScales = options.numScales;
        fprintf('numScales = %d.\n', numScales);
    end
    if isfield(options,'blockSize') && ~isempty(options.blockSize) && isscalar(options.blockSize) && isnumeric(options.blockSize) && options.blockSize > 0
        blockSize = options.blockSize;
        fprintf('blockSize = %d.\n', blockSize);
    end
    if isfield(options,'blockStep') && ~isempty(options.blockStep) && isscalar(options.blockStep) && isnumeric(options.blockStep) && options.blockStep > 0
        blockStep = options.blockStep;
        fprintf('blockStep = %d.\n', blockStep);
    end
    if isfield(options,'hsvBins') && ~isempty(options.hsvBins) && isvector(options.blockStep) && isnumeric(options.hsvBins) && length(options.hsvBins) == 3 && all(options.hsvBins > 0)
        hsvBins = options.hsvBins;
        fprintf('hsvBins = [%d, %d, %d].\n', hsvBins);
    end
    if isfield(options,'tau') && ~isempty(options.tau) && isscalar(options.tau) && isnumeric(options.tau) && options.tau > 0
        tau = options.tau;
        fprintf('tau = %g.\n', tau);
    end
    if isfield(options,'R') && ~isempty(options.R) && isnumeric(options.R) && all(options.R > 0)
        R = options.R;
        fprintf('R = %d.\n', R);
    end
    if isfield(options,'numPoints') && ~isempty(options.numPoints) && isscalar(options.numPoints) && isnumeric(options.numPoints) && options.numPoints > 0
        numPoints = options.numPoints;
        fprintf('numPoints = %d.\n', numPoints);
    end
end


t0 = tic;

timestart = round(timestart * (size(images,4)-1) + 1);
timeend = round(timeend * (size(images,4)-1) + 1);
keyall = round(keynum * (size(images,4)-1) + 1);


descriptors_silpt = [];
descriptors_hsv = [];
descriptors_silpt_video = [];
descriptors_hsv_video = [];
for indscale = 1 : numScales
    feat_hsv = [];
    feat_silpt = [];
    
    for i = 1 : size(images,4)        
        % superpixel
        image = images(:,:,:,i);
        mask = masks(:,:,i);

        % color transformation
        label = SLIC_mex(image, myopt.SLIC_sp_num*3, myopt.SLIC_spatial_proximity_weight);
        N_superpixels = unique(label);
        sp_num = max(N_superpixels); 
        
        % fusion of superpixel and person segmentation
        spmask = [];
        for indNum = 1:sp_num            
            temp_pixel_ratio_mask = sum(mask(find(label==indNum))>0) / length(find(label==indNum));
            spmask(1,indNum) = temp_pixel_ratio_mask > 0.8; 
            [rows cols] = find(label==indNum);
            spmask(2,indNum) = mean(rows);
            spmask(3,indNum) = mean(cols);        
        end

        % extract Joint HSV based LOMO descriptors
        I = Retinex(image);
        imhsv = rgb2hsv(I); 
        feat_hsv = [feat_hsv PyramidMaxJointHist_sp_mask( imhsv, label, spmask, blockSize, blockStep, hsvBins )];       
        feat_silpt = [feat_silpt PyramidMaxSILTPHist_sp_mask( image, label, spmask, blockSize, blockStep, tau, R, numPoints )];
    end
    
    descriptors_hsv_video = [descriptors_hsv_video; feat_hsv];
    descriptors_silpt_video = [descriptors_silpt_video; feat_silpt];
    
    % next scale
    if indscale < numScales
        images = uint8(ColorPooling(images, 'average'));
        mk = double(masks>0);
        mk = Pooling(mk, 'average');
        masks = mk>=0.5;
    end
end

if strcmp(poolingmethod,'avg')
    for i = 1:size(timestart,2)    
        descriptors_hsv = [descriptors_hsv ; mean(descriptors_hsv_video(:,timestart(i):timeend(i)),2)];         %mean
        descriptors_silpt = [descriptors_silpt ; mean(descriptors_silpt_video(:,timestart(i):timeend(i)),2)];         %mean
    end
elseif strcmp(poolingmethod,'max')
    for i = 1:size(timestart,2)    
        descriptors_hsv = [descriptors_hsv ; max(descriptors_hsv_video(:,timestart(i):timeend(i)),[],2)];     %max 
        descriptors_silpt = [descriptors_silpt ; max(descriptors_silpt_video(:,timestart(i):timeend(i)),[],2)];     %max         
    end
elseif strcmp(poolingmethod,'key')
    for i = 1:size(keyall,2)    
        descriptors_hsv = [descriptors_hsv ; descriptors_hsv_video(:,keyall(i))];        %key
        descriptors_silpt = [descriptors_silpt ; descriptors_silpt_video(:,keyall(i))];        %key
    end
end

descriptors_hsv = log(descriptors_hsv + 1);
descriptors_hsv = normc(descriptors_hsv);

descriptors_silpt = log(descriptors_silpt + 1);
descriptors_silpt = normc(descriptors_silpt);

descriptors = [descriptors_hsv; descriptors_silpt];

feaTime = toc(t0);
meanTime = feaTime / size(images, 4);
fprintf('LOMO feature extraction finished. Running time: %.3f seconds in total, %.3f seconds per image.\n', feaTime, meanTime);
end


function descriptors = PyramidMaxJointHist_sp_mask( image, label, spmask, blockSize, blockStep, colorBins )
%% PyramidMaxJointHist: HSV based LOMO representation
    if nargin == 3
        blockSize = 10;
        blockStep = 5;
        colorBins = [8,8,8];
    end

    totalBins = prod(colorBins);  
    
    I = image;

    I(:,:,1) = min( floor( I(:,:,1) * colorBins(1) ), colorBins(1)-1 );
    I(:,:,2) = min( floor( I(:,:,2) * colorBins(2) ), colorBins(2)-1 );
    I(:,:,3) = min( floor( I(:,:,3) * colorBins(3) ), colorBins(3)-1 );
    I = uint16(I);
 
    descriptors = [];
    
    % color feature
    pattern = I(:,:,3) * colorBins(2) * colorBins(1) + I(:,:,2)*colorBins(1) + I(:,:,1); % HSV
    
    % Scan multi-scale blocks and compute histograms
    height = size(image, 1);
    minRow = 1;    
    maxRow = height - blockSize + 1;
    striprow = minRow:blockStep:maxRow;
    stipenum = length(striprow); 
        
    for stnum = 1:stipenum
        temp_feat_strip = [];
        spnum = find(spmask(1,:)==1 & (spmask(2,:) >= striprow(stnum) & spmask(2,:) < striprow(stnum)+blockSize ) );
%         spnum = find( (spmask(2,:) >= striprow(stnum) & spmask(2,:) < striprow(stnum)+blockSize ) );

        for k = 1 : length(spnum)
            temp_sp = pattern(find(label==spnum(k)));       
%             temp_sp = hist(temp_sp, 0 : totalBins-1) ./ length(temp_sp);
            temp_sp = hist(temp_sp, 0 : totalBins-1) ;
            temp_feat_strip = [temp_feat_strip temp_sp'];
        end
        
        if size(temp_feat_strip,2)>0
            descriptors = [descriptors; max(temp_feat_strip,[],2)];
        else
            descriptors = [descriptors; zeros(totalBins,1)];
        end
    end    
end


function descriptors = PyramidMaxJointHist_patch_mask( oriImgs, masks, numScales, blockSize, blockStep, colorBins )
%% PyramidMaxJointHist: HSV based LOMO representation
    if nargin == 2
        numScales = 3;
        blockSize = 10;
        blockStep = 5;
        colorBins = [8,8,8];
    end

    totalBins = prod(colorBins);
    numImgs = size(oriImgs, 4);
    images = zeros(size(oriImgs));

    % color transformation
    for i = 1 : numImgs
        I = oriImgs(:,:,:,i);
        I = Retinex(I);

        I = rgb2hsv(I);
        I(:,:,1) = min( floor( I(:,:,1) * colorBins(1) ), colorBins(1)-1 );
        I(:,:,2) = min( floor( I(:,:,2) * colorBins(2) ), colorBins(2)-1 );
        I(:,:,3) = min( floor( I(:,:,3) * colorBins(3) ), colorBins(3)-1 );
        images(:,:,:,i) = I; % HSV
    end

    minRow = 1;
    minCol = 1;
    descriptors = [];

    % Scan multi-scale blocks and compute histograms
    for i = 1 : numScales
        patterns = images(:,:,3,:) * colorBins(2) * colorBins(1) + images(:,:,2,:)*colorBins(1) + images(:,:,1,:); % HSV
        % changxin
        patterns = squeeze(patterns);
        patterns = patterns .* double(masks>0);
        patterns = reshape(patterns, [], numImgs);
        
        % changxin
        patterns_mask = masks(:,:,:)>0;
        patterns_mask = reshape(patterns_mask, [], numImgs);
        
        height = size(images, 1);
        width = size(images, 2);
        maxRow = height - blockSize + 1;
        maxCol = width - blockSize + 1;

        [cols,rows] = meshgrid(minCol:blockStep:maxCol, minRow:blockStep:maxRow); % top-left positions
        cols = cols(:);
        rows = rows(:);
        numBlocks = length(cols);
        numBlocksCol = length(minCol:blockStep:maxCol);

        if numBlocks == 0
            break;
        end

        offset = bsxfun(@plus, (0 : blockSize-1)', (0 : blockSize-1) * height); % offset to the top-left positions. blockSize-by-blockSize
        index = sub2ind([height, width], rows, cols);
        index = bsxfun(@plus, offset(:), index'); % (blockSize*blockSize)-by-numBlocks
        patches = patterns(index(:), :); % (blockSize * blockSize * numBlocks)-by-numImgs
        patches = reshape(patches, [], numBlocks * numImgs); % (blockSize * blockSize)-by-(numBlocks * numChannels * numImgs)
        fea = hist(patches, 0 : totalBins-1); % totalBins-by-(numBlocks * numImgs)

        % changxin
        patches_mask = patterns_mask(index(:), :);
        patches_mask = reshape(patches_mask, [], numBlocks * numImgs);
        fea_mask = hist(double(patches_mask),0:1);

        fea(1,:) = fea(1,:) - fea_mask(1,:) ;
%         fea(1,:) = 0 ;
%         fea(:,find(fea_mask(1,:)>blockSize*blockSize*0.85)) = 0;
        
        fea = reshape(fea, [totalBins, numBlocks / numBlocksCol, numBlocksCol, numImgs]);
        fea = max(fea, [], 3);
        fea = reshape(fea, [], numImgs);

        descriptors = [descriptors; fea]; %#ok<AGROW>

        if i < numScales
            images = ColorPooling(images, 'average');
            % changxin
            mk = double(masks>0);
            mk = Pooling(mk, 'average');
            masks = mk>=0.5;
%             masks = Pooling(masks, 'average');
        end
    end

%     descriptors = log(descriptors + 1);
%     descriptors = normc(descriptors); 
end


function outImages = ColorPooling(images, method)
    [height, width, numChannels, numImgs] = size(images);
    outImages = images;
    
    if mod(height, 2) == 1
        outImages(end, :, :, :) = [];
        height = height - 1;
    end
    
    if mod(width, 2) == 1
        outImages(:, end, :, :) = [];
        width = width - 1;
    end
    
    if height == 0 || width == 0
        error('Over scaled image: height=%d, width=%d.', height, width);
    end
    
    height = height / 2;
    width = width / 2;
    
    outImages = reshape(outImages, 2, height, 2, width, numChannels, numImgs);
    outImages = permute(outImages, [2, 4, 5, 6, 1, 3]);
    outImages = reshape(outImages, height, width, numChannels, numImgs, 2*2);
    
    if strcmp(method, 'average')
        outImages = floor(mean(outImages, 5));
    else if strcmp(method, 'max')
            outImages = max(outImages, [], 5);
        else
            error('Error pooling method: %s.', method);
        end
    end
end

function descriptors = PyramidMaxSILTPHist_sp_mask( image, label, spmask, blockSize, blockStep, tau, R, numPoints )
%% PyramidMaxSILTPHist: SILTP based LOMO representation

    if nargin == 3
        blockSize = 10;
        blockStep = 5;
        tau = 0.3;
        R = 5;
        numPoints = 4;
    end

    totalBins = 3^numPoints;

    % Convert gray images  
    I = rgb2gray(image);
    image = double(I) / 255;
        
    descriptors = [];
    
    height = size(image, 1);
    minRow = 1;    
    maxRow = height - blockSize + 1;
    striprow = minRow:blockStep:maxRow;
    stipenum = length(striprow); 

    
        
    for stnum = 1:stipenum        
        spnum = find(spmask(1,:)==1 & (spmask(2,:) >= striprow(stnum) & spmask(2,:) < striprow(stnum)+blockSize ) );
%         spnum = find( (spmask(2,:) >= striprow(stnum) & spmask(2,:) < striprow(stnum)+blockSize ) );
        for j = 1: length(R)
            temp_feat_strip = [];
            % siltp feature
            pattern = SILTP(image, tau, R(j), numPoints);
            for k = 1 : length(spnum)
                temp_sp = pattern(find(label==spnum(k)));       
                %             temp_sp = hist(temp_sp, 0 : totalBins-1) ./ length(temp_sp);
                temp_sp = hist(temp_sp, 0 : totalBins-1) ;
                temp_feat_strip = [temp_feat_strip temp_sp'];
            end
        
            if size(temp_feat_strip,2)>0
                descriptors = [descriptors; max(temp_feat_strip,[],2)];
            else
                descriptors = [descriptors; zeros(totalBins,1)];
            end
            
        end
    end     
end

function descriptors = PyramidMaxSILTPHist_patch_mask( oriImgs, masks, numScales, blockSize, blockStep, tau, R, numPoints )
%% PyramidMaxSILTPHist: SILTP based LOMO representation

    if nargin == 2
        numScales = 3;
        blockSize = 10;
        blockStep = 5;
        tau = 0.3;
        R = 5;
        numPoints = 4;
    end

    totalBins = 3^numPoints;

    [imgHeight, imgWidth, ~, numImgs] = size(oriImgs);
    images = zeros(imgHeight,imgWidth, numImgs);

    % Convert gray images
    for i = 1 : numImgs
        I = oriImgs(:,:,:,i);
        I = rgb2gray(I);
        images(:,:,i) = double(I) / 255;
        
    end

    minRow = 1;
    minCol = 1;
    descriptors = [];

    % Scan multi-scale blocks and compute histograms
    for i = 1 : numScales
        height = size(images, 1);
        width = size(images, 2);
        
        if width < R * 2 + 1
            fprintf('Skip scale R = %d, width = %d.\n', R, width);
            continue;
        end
        
        patterns = SILTP(images, tau, R, numPoints);
        % changxin
        patterns = patterns .* double(masks>0);
        patterns = reshape(patterns, [], numImgs);
        
        % changxin
        patterns_mask = masks(:,:,:)>0;
        patterns_mask = reshape(patterns_mask, [], numImgs);
        
        maxRow = height - blockSize + 1;
        maxCol = width - blockSize + 1;

        [cols,rows] = meshgrid(minCol:blockStep:maxCol, minRow:blockStep:maxRow); % top-left positions
        cols = cols(:);
        rows = rows(:);
        numBlocks = length(cols);
        numBlocksCol = length(minCol:blockStep:maxCol);

        if numBlocks == 0
            break;
        end

        offset = bsxfun(@plus, (0 : blockSize-1)', (0 : blockSize-1) * height); % offset to the top-left positions. blockSize-by-blockSize
        index = sub2ind([height, width], rows, cols);
        index = bsxfun(@plus, offset(:), index'); % (blockSize*blockSize)-by-numBlocks
        patches = patterns(index(:), :); % (blockSize * blockSize * numBlocks)-by-numImgs
        patches = reshape(patches, [], numBlocks * numImgs); % (blockSize * blockSize)-by-(numBlocks * numChannels * numImgs)
        
        fea = hist(patches, 0:totalBins-1); % totalBins-by-(numBlocks * numImgs)
        
        % changxin
        patches_mask = patterns_mask(index(:), :);
        patches_mask = reshape(patches_mask, [], numBlocks * numImgs);
        fea_mask = hist(double(patches_mask),0:1);

        fea(1,:) = fea(1,:) - fea_mask(1,:) ;
%         fea(1,:) = 0 ;
%         fea(:,find(fea_mask(1,:)>blockSize*blockSize*0.85)) = 0;

        fea = reshape(fea, [totalBins, numBlocks / numBlocksCol, numBlocksCol, numImgs]);
        fea = max(fea, [], 3);
        fea = reshape(fea, [], numImgs);
        
        descriptors = [descriptors; fea]; %#ok<AGROW>

        if i < numScales
            images = Pooling(images, 'average');
            % changxin
            mk = double(masks>0);
            mk = Pooling(mk, 'average');
            masks = mk>=0.5;
%             masks = Pooling(masks, 'average');
        end
    end

    descriptors = log(descriptors + 1);
    descriptors = normc(descriptors);
end


function outImages = Pooling(images, method)
    [height, width, numImgs] = size(images);
    outImages = images;
    
    if mod(height, 2) == 1
        outImages(end, :, :) = [];
        height = height - 1;
    end
    
    if mod(width, 2) == 1
        outImages(:, end, :) = [];
        width = width - 1;
    end
    
    if height == 0 || width == 0
        error('Over scaled image: height=%d, width=%d.', height, width);
    end
    
    height = height / 2;
    width = width / 2;
    
    outImages = reshape(outImages, 2, height, 2, width, numImgs);
    outImages = permute(outImages, [2, 4, 5, 1, 3]);
    outImages = reshape(outImages, height, width, numImgs, 2*2);
    
    if strcmp(method, 'average')
        outImages = mean(outImages, 4);
    else if strcmp(method, 'max')
            outImages = max(outImages, [], 4);
        else
            error('Error pooling method: %s.', method);
        end
    end
end
