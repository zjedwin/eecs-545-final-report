function cluster_with_GMM(X,k)
% Clusters the data matrix X using GMM after applying PCA
% 'X'   - data matrix with rows as data points
% 'k'   - number of clusters
% Zamar Edwin, Charles Lu

% Y is the representation of X in the PCA space
[~,Y] = pca(X,'NumComponents',3);

% GMM Clustering
gm = fitgmdist(Y,k);
idx = cluster(gm,Y);
gmClusters = zeros(size(Y));

for i=1:k
    % gmClusters(:,n*(i-1)+1:n*i) contains the elements of Y
    % belonging to the i-th cluster.
    n = sum(idx == i);
    gmClusters(n*(i-1)+1:n*i,:) = Y(idx == i,:);
end

% calculate cluster centers
gmCenters = zeros(k,size(X,2));
for c = 1:k
    gmCenters(c,:) = mean(X(idx==c,:));
end

% saves cluster centers and labels to the workspace
assignin('caller','gmmCenters',gmCenters);
assignin('caller','gmmClusterLabels',idx);
end