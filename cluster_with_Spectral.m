function cluster_with_Spectral(X,k,sigma)
% Clusters the data matrix X using spectral clustering after applying PCA
% 'X'       - data matrix with rows as data points
% 'k'       - number of clusters
% 'sigma'   - parameter for Gaussian function in deciding the graph weights

% Y is the representation of X in the PCA space
%[~,Y] = pca(X,'NumComponents',2);

% Zamar Edwin; Charles Lu

A = SimGraph_Full(X',sigma);
[spectralClusterLabels,~,~] = SpectralClustering(A,k,2); % 2 == Normalized
spectralCenters = zeros(k,size(X,2));
for c = 1:k
spectralCenters(c,:) = mean(X(spectralClusterLabels==c,:),1);
end

% saves cluster centers and labels to the workspace
assignin('caller','spectralCenters',spectralCenters);
assignin('caller','spectralClusterLabels',spectralClusterLabels);

end
