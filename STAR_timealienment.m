function [walking_cycle cycle_index]= STAR_timealienment(myopt,cam_num)
%% function STAR_timealienment( myopt, cam_num, cycle_flag )
% Function for walking cycle selection, 
% e.g., time alienment by motion information (position (x) curves) processing
%
% Input:
%   <myopt>: parameters
%   <myopt>: index of camera view
%
% Output:
%   <video_clip>: selected walking cycle
%
% Date: 2016-07-18
% Author: Changxin Gao
% Institute: School of Automation, Huazhong University of Science and Technology
% Email: cgao@hust.edu.cn

%% parameters
camName = myopt.camName{cam_num};
cycle_flag = myopt.cycle_flag;
% load positions
load(['positions\pos_' camName(1:end-1) '.mat']);

walking_cycle = [];
cycle_index = [];
num_ind = 1;

%%
for i = 1:length(person)
    disp(i)
    close all
    posDiff = person(i).posDiff;
    x = 1:size(posDiff,2);
    y_up = x*0 + myopt.col/2 + myopt.threshold;
    y_down = x*0 + myopt.col/2 - myopt.threshold;
        
    clip = [];
    scores = [];
    for j = 1:size(posDiff,1)
        y = posDiff(j,:); 
        
        % curve fitting
        c = polyfit(x, y, 20);  
        d = polyval(c, x, 1);        
        
        if cycle_flag == 0.5
            [score locs pks] = halfCycle(d, x, myopt);  
            if length(locs) > 0
                clip = [clip; reshape(locs,[length(locs)/2 2])];
                scores = [scores; score];
            end        
        elseif cycle_flag == 1
            [score locs pks] = oneCycle(d, x, myopt); 
            if length(locs) > 0
                clip = [clip; reshape(locs,[length(locs)/3 3])];
                scores = [scores; score];
            end
        else
            reture;    
        end
    end

    if strcmp(myopt.numcycle,'multiple')
        s_index = find(scores>0.5);
        if length(s_index)>0        
            for iii = 1:length(s_index)
               walking_cycle(num_ind,:) = clip(iii,:);
               cycle_index(num_ind,:) = i;
               num_ind = num_ind + 1;
            end        
        else if cycle_flag == 1
                walking_cycle(num_ind,:) = [1 15 30];
                cycle_index(num_ind,:) = i;
                num_ind = num_ind + 1;
            elseif cycle_flag == 0.5
                walking_cycle(num_ind,:) = [1 15];
                cycle_index(num_ind,:) = i;
                num_ind = num_ind + 1;
            end
        end  
    else %if myopt.numcycle == 'singl'
        cycle_index(i,:) = i;
        if size(clip,1)>0
            iii = find(scores==max(scores));
            if length(iii)>1
                iii = iii(1);
            end
            walking_cycle(i,:) = clip(iii,:);
        else if cycle_flag == 1
                walking_cycle(i,:) = [1 15 30];
            elseif cycle_flag == 0.5
                walking_cycle(i,:) = [1 15];
            end
        end
    end
    
end
end

function [score locs pks] = halfCycle(d, x, myopt)
Lmax = diff(sign(diff(d)))== -2; % logic vector for the local max value
Lmin = diff(sign(diff(d)))== 2; % logic vector for the local min value
% match the logic vector to the original vecor to have the same length
Lmax = [false Lmax false];
Lmin =  [false Lmin false];
locs_max = x(Lmax); % locations of the local max elements
locs_min = x(Lmin); % locations of the local min elements
pks_max = d(Lmax); % values of the local max elements
pks_min = d(Lmin); % values of the local min elements

pks = [pks_max pks_min];
locs = [locs_max locs_min];
flag = zeros(1,length(pks))-1;
flag(1:length(pks_max))=1;% flag of max (1) or min (0)
[locs_sort ind_sort] = sort(locs);
pks_sort=pks(ind_sort);
flag_sort = flag(ind_sort);   


pks = [];
locs = [];
score = [];
for i = 1:length(pks_sort)-2
    if flag_sort(i)*pks_sort(i) > flag_sort(i)*(myopt.col/2 + flag_sort(i)*myopt.threshold)
        if flag_sort(i+1)*pks_sort(i+1) > flag_sort(i+1)*(myopt.col/2 + flag_sort(i+1)*myopt.threshold)
            locs = [locs locs_sort(i) locs_sort(i+1)];
            pks = [pks pks_sort(i) pks_sort(i+1)]; 
            
            clip = [locs_sort(i) locs_sort(i+1)];               
            
            
            func_sin = @(a,t) a(1)*sin(2*pi*(t-clip(1))/(clip(2)-clip(1))/2 - 0.5*pi) + a(2);                
            A= lsqcurvefit( func_sin, [1 myopt.col/2], clip(1):clip(2), d(clip(1):clip(2)));
         
            temp = locs_sort(i+1) - locs_sort(i);
            if temp > 9 & temp < 22
                score_cycle = 1;
            else
                score_cycle = 0;
            end
            
            weight = d(clip(1):clip(2)) - func_sin(A,clip(1):clip(2));
            weight = 1 - sqrt(sum(weight.^2)/(clip(2)-clip(1)))/myopt.col;
            score = [score; score_cycle + log(weight)];
        end
    end
end

end




function [score locs pks] = oneCycle(d, x, myopt)
Lmax = diff(sign(diff(d)))== -2; % logic vector for the local max value
Lmin = diff(sign(diff(d)))== 2; % logic vector for the local min value
% match the logic vector to the original vecor to have the same length
Lmax = [false Lmax false];
Lmin =  [false Lmin false];
locs_max = x(Lmax); % locations of the local max elements
locs_min = x(Lmin); % locations of the local min elements
pks_max = d(Lmax); % values of the local max elements
pks_min = d(Lmin); % values of the local min elements

pks = [pks_max pks_min];
locs = [locs_max locs_min];
flag = zeros(1,length(pks))-1;
flag(1:length(pks_max))=1;% flag of max (1) or min (0)
[locs_sort ind_sort] = sort(locs);
pks_sort=pks(ind_sort);
flag_sort = flag(ind_sort);   


pks = [];
locs = [];
score = [];    

for i = 1:length(pks_sort)-2
    if flag_sort(i)*pks_sort(i) > flag_sort(i)*(myopt.col/2 + flag_sort(i)*myopt.threshold)
        if flag_sort(i+1)*pks_sort(i+1) > flag_sort(i+1)*(myopt.col/2 + flag_sort(i+1)*myopt.threshold)
            if flag_sort(i+2)*pks_sort(i+2) > flag_sort(i+2)*(myopt.col/2 + flag_sort(i+2)*myopt.threshold)
                locs = [locs locs_sort(i) locs_sort(i+1) locs_sort(i+2)];
                pks = [pks pks_sort(i) pks_sort(i+1) pks_sort(i+2)]; 
                
                clip = [locs_sort(i) locs_sort(i+1) locs_sort(i+2)];
                
                func_sin = @(a,t) a(1)*sin(2*pi*(t-clip(1))/(clip(3)-clip(1)) - 0.5*pi) + a(2);                
                A= lsqcurvefit( func_sin, [1 myopt.col/2], clip(1):clip(3), d(clip(1):clip(3))); 
                
                temp = locs_sort(i+2) - locs_sort(i);
                if temp > 18 & temp < 42
                    score_cycle = 1;
                else
                    score_cycle = 0;
                end
                
                weight = d(clip(1):clip(3)) - func_sin(A,clip(1):clip(3));
                weight = 1 - sqrt(sum(weight.^2)/(clip(3)-clip(1)))/myopt.col;
                score = [score; score_cycle*(1+log(weight))];
            end
        end
    end
end
end