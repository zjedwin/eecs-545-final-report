function cluster_with_SSC(X,k,lambda)
% Clusters the data matrix X using sparse subspace clustering
% 'X'       - data matrix with rows as data points
% 'k'       - number of clusters
% 'sigma'   - parameter for Gaussian function in deciding the graph weights
% Zamar Edwin; Charles Lu

[C,~] = ALM_noisySSC(X',lambda);
W = abs(C) + abs(C');
[sscClusterLabels,~,~] = SpectralClustering(W,k,2); % 2 == Normalized
sscCenters = zeros(k,size(X,2));
for c = 1:k
sscCenters(c,:) = mean(X(sscClusterLabels==c,:),1);
end

% saves cluster centers and labels to the workspace
assignin('caller','sscCenters',sscCenters);
assignin('caller','sscClusterLabels',sscClusterLabels);

end