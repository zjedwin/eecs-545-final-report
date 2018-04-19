function cluster_with_SC(X,k)
% Clusters the data matrix X using spectral clustering after applying PCA
% 'X'       - data matrix with rows as data points
% 'k'       - number of clusters
% 'sigma'   - parameter for Gaussian function in deciding the graph weights
% Charles Lu; Zamar Edwin

% Unweighted adjacency graph
W = SimGraph_NearestNeighbors(X',k,1,0); % 1 == Normal KNN
[SC_labels,SC_centers] = SpectralClustering(W,k,1); % 1 == Unnormalized

[i,j,~] = find(SC_labels);
SC_labels = [i,j];
SC_labels = sortrows(SC_labels);
SC_labels = SC_labels(:,2);

% saves cluster centers and labels to the workspace
assignin('caller','SC_centers',SC_centers);
assignin('caller','SC_labels',SC_labels);

end
