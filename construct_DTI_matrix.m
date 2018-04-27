function M = construct_DTI_matrix(net)
% convert bipartite network to adjacency matrix
m = net(:,1); % drug
d = net(:,2); % target

M = zeros(max(m),max(d));

for i = 1 : length(m)
    M(m(i),d(i)) = 1;
end







