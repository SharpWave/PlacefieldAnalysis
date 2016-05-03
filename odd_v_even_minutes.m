function [odd_ind, even_ind ] = odd_v_even_minutes(frames, SR)
% [odd_ind, even_ind ] = odd_v_even_minutes(frames, SR)
%   Takes a vector frames obtained at sampling rate SR (in frames/sec) and
%   spits out vectors of indices to the odd-minute frames and the even
%   minute frames


duration = length(frames)/SR; % in seconds
time_adj = (1/SR:1/SR:length(frames)/SR)/60; % in minutes

% Generate odd-even limits
odd_lims = repmat(0:2:duration,2,1) + [zeros(1,floor(duration/2)+1); ones(1,floor(duration/2)+1)];
even_lims = repmat(1:2:duration+1,2,1) + [zeros(1,floor(duration/2)+1); ones(1,floor(duration/2)+1)];

odd_ind = false(size(frames));
even_ind = false(size(frames));

%%
for j = 1:size(odd_lims,2)
   odd_ind = odd_ind | (time_adj >= odd_lims(1,j) & time_adj < odd_lims(2,j));
   even_ind = even_ind | (time_adj >= even_lims(1,j) & time_adj < even_lims(2,j));
    
end

odd_ind = find(odd_ind);
even_ind = find(even_ind);

end

