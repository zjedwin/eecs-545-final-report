function cluster_with_kMeans(X,k)
% Clusters the data matrix X using kMeans after applying PCA
% 'X'   - data matrix with rows as data points
% 'k'   - number of clusters
% Zamar Edwin; Charles Lu

% Y is the representation of X in the PCA space
[~,Y] = pca(X,'NumComponents',3);

% kMeans Clustering
[kmClusterLabels,~] = kmeans(Y,k);

% calculate cluster centers
kmCenters = zeros(k,size(X,2));
for c = 1:k
    kmCenters(c,:) = mean(X(kmClusterLabels==c,:));
end

% saves cluster centers and labels to the workspace
assignin('caller', 'kmCenters', kmCenters);
assignin('caller', 'kmClusterLabels', kmClusterLabels);
end