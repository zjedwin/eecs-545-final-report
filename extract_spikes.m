%% EECS 545 project
% extracts spikes from simulated neural recording at
% http://www2.le.ac.uk/departments/engineering/research/bioengineering/neuroengineering-lab/software
% Zamar Edwin; Charles Lu

%% Parameters
T_s = 64;   % number of samples to extract per spike
T_p = 30;   % number of samples preceding (and including) spike peak

%% Process data and parameters
i_spikes = spike_times{1};
N = length(i_spikes);

%% Extract spikes
X = zeros(N,T_s);
for n = 1:N
    snippet = data(i_spikes(n):i_spikes(n)+T_s);
    [~,i_max] = max(snippet);
    i_center = i_spikes(n) + i_max - 1;
    X(n,:) = data(i_center-T_p:i_center+T_s-T_p-1);
end
