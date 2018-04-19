%% EECS 545 project
% Philip Vu

% GMM SimGaussian data
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\gm_gmm_outputs.mat')
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Ready-to-use data sets\DATA-mixture_gaussians.mat')

indx = 0;
MaxindxPerCluster = zeros(1,5);
indxsPerClus = cell(1,5);
for i = [-1 1 2 3 4]
    indx = indx + 1; 
    indxsPerClus{indx} = find(Y == i); 
    MaxindxPerCluster(indx) = max(indxsPerClus{indx});
end

[trash tempIndx] = sort([MaxindxPerCluster], 'descend');
indxsPerClus = indxsPerClus(tempIndx);
%% 
MaxindxPerClusterGMM = zeros(10,5);
indxsPerClusGMM = cell(10,5); 
for i = 1:length(gm_gmm_outputs)
    for j = 1:5
        indxsPerClusGMM{i,j} = find(gm_gmm_outputs(i).labels == j);
        MaxindxPerClusterGMM(i,j) = max(indxsPerClusGMM{i,j});
    end
    [trash tempIndx] = sort([MaxindxPerClusterGMM(i,:)], 'descend');
    indxsPerClusGMM(i,:) = indxsPerClusGMM(i,tempIndx);
end
%% 

clusAcc = zeros(10,5);
nMissClassified = zeros(10,5); 
for i = 1:10 
    for j = 1:5
        clusAcc(i,j) = (length(indxsPerClusGMM{i,j})-length(setdiff(indxsPerClus{1,j},indxsPerClusGMM{i,j})))/length(indxsPerClusGMM{i,j}); %Calculate accuracy of cluster algorithm
        nMissClassified(i,j) = length(setdiff(indxsPerClusGMM{i,j},indxsPerClus{1,j})); %Count number of misclassified datapoints 
    end
end 
avgClusAcc = mean(clusAcc)
avgMissClassified = mean(nMissClassified) 


%% K-Means simgaussian data
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\gm_km_outputs.mat')
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Ready-to-use data sets\DATA-mixture_gaussians.mat')

indx = 0;
MaxindxPerCluster = zeros(1,5);
indxsPerClus = cell(1,5);
for i = [-1 1 2 3 4]
    indx = indx + 1;
    indxsPerClus{indx} = find(Y == i);
    MaxindxPerCluster(indx) = max(indxsPerClus{indx});
end

[trash tempIndx] = sort([MaxindxPerCluster], 'descend');
indxsPerClus = indxsPerClus(tempIndx);
%% 
MaxindxPerClusterKM = zeros(10,5);
indxsPerClusKM = cell(10,5); 
for i = 1:length(gm_km_outputs)
    for j = 1:5
        indxsPerClusKM{i,j} = find(gm_km_outputs(i).labels == j);
        MaxindxPerClusterKM(i,j) = max(indxsPerClusKM{i,j});
    end
    [trash tempIndx] = sort([MaxindxPerClusterKM(i,:)], 'descend');
    indxsPerClusKM(i,:) = indxsPerClusKM(i,tempIndx);
end
tempIndx = [1 4 5 3 2];
indxsPerClusKM(i,:) = indxsPerClusKM(i,tempIndx);
%% 

clusAcc = zeros(10,5);
nMissClassified = zeros(10,5); 
for i = 1:10 
    for j = 1:5
        clusAcc(i,j) = (length(indxsPerClusKM{i,j})-length(setdiff(indxsPerClus{1,j},indxsPerClusKM{i,j})))/length(indxsPerClusKM{i,j}); %Calculate accuracy of cluster algorithm
        nMissClassified(i,j) = length(setdiff(indxsPerClusKM{i,j},indxsPerClus{1,j})); %Count number of misclassified datapoints 
    end
end 
avgClusAcc = mean(clusAcc)
avgMissClassified = mean(nMissClassified) 
%% Spectral simgaussian data
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\gm_spectral_outputs.mat')
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Ready-to-use data sets\DATA-mixture_gaussians.mat')

indx = 0;
MaxindxPerCluster = zeros(1,5);
indxsPerClus = cell(1,5);
for i = [-1 1 2 3 4]
    indx = indx + 1;
    indxsPerClus{indx} = find(Y == i);
    MaxindxPerCluster(indx) = max(indxsPerClus{indx});
end

[trash tempIndx] = sort([MaxindxPerCluster], 'descend');
indxsPerClus = indxsPerClus(tempIndx);
%% 
MaxindxPerClusterSC = zeros(1,5);
indxsPerClusSC = cell(1,5); 
for i = 1
    for j = 1:5
        indxsPerClusSC{i,j} = find(spectralClusterLabels == j);
        MaxindxPerClusterSC(i,j) = max(indxsPerClusSC{i,j});
    end
    [trash tempIndx] = sort([MaxindxPerClusterSC(i,:)], 'descend');
    indxsPerClusSC(i,:) = indxsPerClusSC(i,tempIndx);
end
% tempIndx = [1 4 5 3 2];
% indxsPerClusKM(i,:) = indxsPerClusKM(i,tempIndx);
%% 

clusAcc = zeros(1,5);
nMissClassified = zeros(1,5); 
for i = 1
    for j = 1:5
        clusAcc(i,j) = (length(indxsPerClusSC{i,j})-length(setdiff(indxsPerClus{1,j},indxsPerClusSC{i,j})))/length(indxsPerClusSC{i,j}); %Calculate accuracy of cluster algorithm
        nMissClassified(i,j) = length(setdiff(indxsPerClusSC{i,j},indxsPerClus{1,j})); %Count number of misclassified datapoints 
    end
end 
avgClusAcc = mean(clusAcc)
avgMissClassified = mean(nMissClassified) 
%% SSC simgaussian data
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\gm_ssc_outputs.mat')
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Ready-to-use data sets\DATA-mixture_gaussians.mat')

indx = 0;
MaxindxPerCluster = zeros(1,5);
indxsPerClus = cell(1,5);
for i = [-1 1 2 3 4]
    indx = indx + 1;
    indxsPerClus{indx} = find(Y == i);
    MaxindxPerCluster(indx) = max(indxsPerClus{indx});
end

[trash tempIndx] = sort([MaxindxPerCluster], 'descend');
indxsPerClus = indxsPerClus(tempIndx);
%% 
MaxindxPerClusterSSC = zeros(1,5);
indxsPerClusSSC = cell(1,5); 
for i = 1
    for j = 1:5
        indxsPerClusSSC{i,j} = find(sscClusterLabels == j);
        MaxindxPerClusterSSC(i,j) = max(indxsPerClusSSC{i,j});
    end
    [trash tempIndx] = sort([MaxindxPerClusterSSC(i,:)], 'descend');
    indxsPerClusSSC(i,:) = indxsPerClusSSC(i,tempIndx);
end
% tempIndx = [1 4 5 3 2];
% indxsPerClusKM(i,:) = indxsPerClusKM(i,tempIndx);
%% 

clusAcc = zeros(1,5);
nMissClassified = zeros(1,5); 
for i = 1
    for j = 1:5
        clusAcc(i,j) = (length(indxsPerClusSSC{i,j})-length(setdiff(indxsPerClus{1,j},indxsPerClusSSC{i,j})))/length(indxsPerClusSSC{i,j}); %Calculate accuracy of cluster algorithm
        nMissClassified(i,j) = length(setdiff(indxsPerClusSSC{i,j},indxsPerClus{1,j})); %Count number of misclassified datapoints 
    end
end 
avgClusAcc = mean(clusAcc)
avgMissClassified = mean(nMissClassified) 
%% Plot GMM Real Spikes
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\realspike_gmm_outputs.mat')

colors = {'b','c','r','g','m','k'};
for j = 2
    figure(j)
    for i = 1:6
        
        indxs = find(gmm_outputs(j).labels == i);
        %subplot(1,4,indx)
        subplot(2,3,i)
        h1 = shadedErrorBar(1:size(X(indxs,:),2), X(indxs,:), {@mean, @std}, {colors{1}, 'LineWidth', 0.1}, 1);
        if(i == 2)
            title('Real Spike Data Clustered Waveforms GMM')
        end
        if(i == 1 || i == 4)
            ylabel('Voltage \muV')
        end
        if(i == 5)
            xlabel('Samples')
        end
        %plot(X(indxs,:)', colors{indx})
        %hold on
    end
%     title('Real Spike Data Clustered Waveforms GMM')
%     ylabel('Voltage \muV')
%     xlabel('Samples')
    %legend('Unit 1', 'Unit 2', 'Unit 3', 'Unit 4', 'Unit 5', 'Unit 6')
end
%% Plot K-Means Real Spikes
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\realspike_km_outputs.mat')
colors = {'b','c','r','g','m','k'};
for j = 9
    figure(j)
    for i = 1:6
        
        indxs = find(km_outputs(j).labels == i);
        %subplot(1,4,indx)
        subplot(2,3,i)
        h1 = shadedErrorBar(1:size(X(indxs,:),2), X(indxs,:), {@mean, @std}, {colors{1}, 'LineWidth', 0.1}, 1);
        if(i == 2)
            title('Real Spike Data Clustered Waveforms K-Means')
        end
        if(i == 1 || i == 4)
            ylabel('Voltage \muV')
        end
        if(i == 5)
            xlabel('Samples')
        end
        %plot(X(indxs,:)', colors{indx})
        %hold on
    end
%     title('Real Spike Data Clustered Waveforms K-Means')
%     ylabel('Voltage \muV')
%     xlabel('Samples')
    %legend('Unit 1', 'Unit 2', 'Unit 3', 'Unit 4', 'Unit 5', 'Unit 6')
end
%% Plot Spectral Real Spikes
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\realspike_spectral_outputs.mat')

colors = {'b','c','r','g','m','k'};
for j = 1
    figure(j)
    for i = 1:6
        
        indxs = find(spectralClusterLabels == i);
        
        subplot(2,3,i)
        h1 = shadedErrorBar(1:size(X(indxs,:),2), X(indxs,:), {@mean, @std}, {colors{1}, 'LineWidth', 0.1}, 1);
        %plot(X(indxs,:)', colors{indx})
        %hold on
        if(i == 2)
            title('Real Spike Data Clustered Waveforms SC')
        end
        if(i == 1 || i == 4)
            ylabel('Voltage \muV')
        end
        if(i == 5)
            xlabel('Samples')
        end
    end
    
end
   
    %legend('Unit 1', 'Unit 2', 'Unit 3', 'Unit 4', 'Unit 5', 'Unit 6')

%% Plot SSC Real Spikes
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\realspike_ssc_outputs.mat')

colors = {'b','c','r','g','m','k'};
for j = 1
    figure(j)
    for i = 1:6
        
        indxs = find(sscClusterLabels == i);
        subplot(2,3,i)
        h1 = shadedErrorBar(1:size(Xd(indxs,:),2), Xd(indxs,:), {@mean, @std}, {colors{1}, 'LineWidth', 0.1}, 1);
        %plot(X(indxs,:)', colors{indx})
        %hold on
        if(i == 2)
            title('Real Spike Data Clustered Waveforms SSC')
        end
        if(i == 1 || i == 4)
            ylabel('Voltage \muV')
        end
        if(i == 5)
            xlabel('Samples')
        end
    end
    
    %legend('Unit 1', 'Unit 2', 'Unit 3', 'Unit 4', 'Unit 5', 'Unit 6')
end
%% GMM analysis
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\simspike_gmm_outputs.mat')
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Ready-to-use data sets\DATA-simulated_spikes.mat')

colors = {'c','r','b','g','m'};
indx = 0;
nSamplesPerCluster = zeros(1,4);
indxsPerClus = cell(1,4);
for i = [-1 0 1 2]
    indx = indx + 1; 
    nSamplesPerCluster(indx) = length(find(Y == i)); 
    indxsPerClus{indx} = find(Y == i); 
    
end
[trash tempIndx] = sort([nSamplesPerCluster], 'descend');
indxsPerClus = indxsPerClus(tempIndx);
%% 
nSamplesPerClusterGMM = zeros(10,4);
indxsPerClusGMM = cell(10,4); 
for i = 1:length(simspike_gmm_outputs)
    for j = 1:4
        nSamplesPerClusterGMM(i,j) = length(find(simspike_gmm_outputs(i).labels == j));
        indxsPerClusGMM{i,j} = find(simspike_gmm_outputs(i).labels == j);
    end
    [trash tempIndx] = sort([nSamplesPerClusterGMM(i,:)], 'descend');
    indxsPerClusGMM(i,:) = indxsPerClusGMM(i,tempIndx);
end
%% 

clusAcc = zeros(10,4);
nMissClassified = zeros(10,4); 
for i = 1:10 
    for j = 1:4
        clusAcc(i,j) = (length(indxsPerClusGMM{i,j})-length(setdiff(indxsPerClus{1,j},indxsPerClusGMM{i,j})))/length(indxsPerClusGMM{i,j});
        nMissClassified(i,j) = length(setdiff(indxsPerClusGMM{i,j},indxsPerClus{1,j}));
    end
end 
avgClusAcc = mean(clusAcc)
avgMissClassified = mean(nMissClassified) 

%% Plot simulated spike data Waveforms GMM 

colors = {'c','r','b','g','m'};
figure(4)
clf
for i = 1:4
    indxs = find(simspike_gmm_outputs(4).labels == i);
    %subplot(1,4,i)
    h1 = shadedErrorBar(1:size(X(indxs,:),2), X(indxs,:), {@mean, @std}, {colors{i}, 'LineWidth', .1}, 1);
    %plot(X(indxs,:)', colors{i})
    hold on
end
title('Simulated Spike Data Clustered Waveforms (GMM)')
ylabel('Voltage \muV')
xlabel('Samples') 

%% K-Means analysis 
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\simspike_km_outputs.mat')
indx = 0;
nSamplesPerCluster = zeros(1,4);
indxsPerClus = cell(1,4);
for i = [-1 0 1 2]
    indx = indx + 1; 
    nSamplesPerCluster(indx) = length(find(Y == i)); 
    indxsPerClus{indx} = find(Y == i); 
   
end
[trash tempIndx] = sort([nSamplesPerCluster], 'descend');
indxsPerClus = indxsPerClus(tempIndx);
%% 
nSamplesPerClusterKM = zeros(10,4);
indxsPerClusKM = cell(10,4); 
for i = 1:length( simspike_km_outputs)
    for j = 1:4
        nSamplesPerClusterKM(i,j) = length(find( simspike_km_outputs(i).labels == j));
        indxsPerClusKM{i,j} = find( simspike_km_outputs(i).labels == j);
    end
    [trash tempIndx] = sort([nSamplesPerClusterKM(i,:)], 'descend');
    indxsPerClusKM(i,:) = indxsPerClusKM(i,tempIndx);
end
%% 

clusAcc = zeros(10,4);
nMissClassified = zeros(10,4); 
for i = 1:10 
    for j = 1:4
        clusAcc(i,j) = (length(indxsPerClusKM{i,j})-length(setdiff(indxsPerClus{1,j},indxsPerClusKM{i,j})))/length(indxsPerClusKM{i,j});
        nMissClassified(i,j) = length(setdiff(indxsPerClusKM{i,j},indxsPerClus{1,j}));
    end
end 
avgClusAcc = mean(clusAcc)
avgMissClassified = mean(nMissClassified) 

%% Plot simulated spike data Waveforms KM 

colors = {'c','r','b','g','m'};
figure(4)
clf
for i = 1:4
    indxs = find(simspike_km_outputs(4).labels == i);
    %subplot(1,4,i)
    h1 = shadedErrorBar(1:size(X(indxs,:),2), X(indxs,:), {@mean, @std}, {colors{i}, 'LineWidth', .1}, 1);
    %plot(X(indxs,:)', colors{i})
    hold on
end
title('Simulated Spike Data Clustered Waveforms (K-Means)')
ylabel('Voltage \muV')
xlabel('Samples') 
%% Spectral Analysis
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\simspike_spectral_outputs.mat')
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Ready-to-use data sets\DATA-simulated_spikes.mat')

indx = 0;
nSamplesPerCluster = zeros(1,4);
indxsPerClus = cell(1,4);
for i = [-1 0 1 2]
    indx = indx + 1; 
    nSamplesPerCluster(indx) = length(find(Y == i)); 
    indxsPerClus{indx} = find(Y == i); 
   
end
[trash tempIndx] = sort([nSamplesPerCluster], 'descend');
indxsPerClus = indxsPerClus(tempIndx);
%% 
nSamplesPerClusterSC = zeros(1,4);
indxsPerClusSC = cell(1,4); 

for j = 1:4
    nSamplesPerClusterSC(j) = length(find( spectralClusterLabels == j));
    indxsPerClusSC{j} = find( spectralClusterLabels == j);
end
[trash tempIndx] = sort([nSamplesPerClusterSC], 'descend');
indxsPerClusSC(1,:) = indxsPerClusSC(1,tempIndx);

%% 

clusAcc = zeros(1,4);
nMissClassified = zeros(1,4); 
for i = 1
    for j = 1:4
        clusAcc(i,j) = (length(indxsPerClusSC{i,j})-length(setdiff(indxsPerClus{1,j},indxsPerClusSC{i,j})))/length(indxsPerClusSC{i,j});
        nMissClassified(i,j) = length(setdiff(indxsPerClusSC{i,j},indxsPerClus{1,j}));
    end
end 
%avgClusAcc = mean(clusAcc)
%avgMissClassified = mean(nMissClassified) 

%% Plot simulated spike data Waveforms Spectral

colors = {'m','c','r','b','c'};
figure(4)
clf
for i = [2 3 4 1]
    indxs = find(spectralClusterLabels == i);
    %subplot(1,4,i)
    h1 = shadedErrorBar(1:size(X(indxs,:),2), X(indxs,:), {@mean, @std}, {colors{i}, 'LineWidth', .1}, 1);
    %plot(X(indxs,:)', colors{i})
    hold on
end
title('Simulated Spike Data Clustered Waveforms (SC)')
ylabel('Voltage \muV')
xlabel('Samples') 
%% SSC Analysis
clear
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Classifier outputs\simspike_ssc_outputs.mat')
load('C:\Users\chesteklab\Documents\2016 Fall\EECS 545 (Machine Learning)\Project\Ready-to-use data sets\DATA-simulated_spikes.mat')

indx = 0;
nSamplesPerCluster = zeros(1,4);
indxsPerClus = cell(1,4);
for i = [-1 0 1 2]
    indx = indx + 1; 
    nSamplesPerCluster(indx) = length(find(Y == i)); 
    indxsPerClus{indx} = find(Y == i); 
   
end
[trash tempIndx] = sort([nSamplesPerCluster], 'descend');
indxsPerClus = indxsPerClus(tempIndx);
%% 
nSamplesPerClusterSSC = zeros(1,4);
indxsPerClusSSC = cell(1,4); 

for j = 1:4
    nSamplesPerClusterSSC(j) = length(find( sscClusterLabels == j));
    indxsPerClusSSC{j} = find( sscClusterLabels == j);
end
[trash tempIndx] = sort([nSamplesPerClusterSSC], 'descend');
indxsPerClusSSC(1,:) = indxsPerClusSSC(1,tempIndx);

%% 

clusAcc = zeros(1,4);
nMissClassified = zeros(1,4); 
for i = 1
    for j = 1:4
        clusAcc(i,j) = (length(indxsPerClusSSC{i,j})-length(setdiff(indxsPerClus{1,j},indxsPerClusSSC{i,j})))/length(indxsPerClusSSC{i,j});
        nMissClassified(i,j) = length(setdiff(indxsPerClusSSC{i,j},indxsPerClus{1,j}));
    end
end 

%% Plot simulated spike data Waveforms SSC

colors = {'c','b','g','m','r'};
figure(4)
clf
for i = [2 3 4 1]
    indxs = find(sscClusterLabels == i);
    %subplot(1,4,i)
    h1 = shadedErrorBar(1:size(X(indxs,:),2), X(indxs,:), {@mean, @std}, {colors{i}, 'LineWidth', .1}, 1);
    %plot(X(indxs,:)', colors{i})
    hold on
end
title('Simulated Spike Data Clustered Waveforms (SSC)')
ylabel('Voltage \muV')
xlabel('Samples') 
%% Plot simulated spike data Original Waveforms
figure(3)
clf
indx = 0; 
colors = {'b','c','r','g','m'};
for i = [-1 0 1 2]
    indx = indx + 1;
    indxs = find(Y == i);
    %subplot(1,4,indx)
    h1 = shadedErrorBar(1:size(X(indxs,:),2), X(indxs,:), {@mean, @std}, {colors{indx}, 'LineWidth', .1}, 1);
   %plot(X(indxs,:)', colors{indx})
    hold on
end
title('Simulated Spike Data Original Waveforms')
ylabel('Voltage \muV')
xlabel('Samples')


