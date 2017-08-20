function k = t1_filename(f, num)
%% Copyright (C) Shu Wang.
%% All rights reserved.
% turn a number smaller than 1,000,000 into a num digit number.
temp_str = int2str(f);
temp_length = length(temp_str);
if(temp_length == num)
    k = int2str(f);
else
    for i = 1:(num - temp_length)
        temp_str = ['0' temp_str];
    end
end
k = [temp_str];
