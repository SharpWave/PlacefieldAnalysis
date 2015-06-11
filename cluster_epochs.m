function [ on_cluster, off_cluster ] =cluster_epochs( on, off, thresh )
% Takes two vectors with indices when a binary vector turns on or off and
% clumps/clusters that are very close to one another (subsequent on indices
% < thresh) and spits them out in on_cluster and off_cluster

n = 1;
temp_on(1) = on(1);
n = 1;
for j = 1:length(on)-1
    if on(j+1) - on(j) < thresh
        continue
    end
    temp_off(n) = off(j);
    n = n + 1;
    temp_on(n) = on(j+1);
end

temp_off(n) = off(end);

on_cluster = temp_on;
off_cluster = temp_off;

end

