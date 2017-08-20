function probability = t1_cal_gaussian_pro(data, cent, sig)
%% Copyright (C) Shu Wang.
%% All rights reserved.
% motion model based on gaussian distribution
dist = slmetric_pw(data, cent, 'sqdist');
probability = (2*pi*sig)^(-.5) * exp(-0.5*dist/sig);
probability = probability';
