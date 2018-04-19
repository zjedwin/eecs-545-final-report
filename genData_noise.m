%% EECS 545 project
% generates pink noise and then picks out thresholded snippets to act as
% noisy signals to make clustering harder
% Charles Lu

%% Parameters
N = 1000;   % number of snippets to create
T_s = 64;   % number of samples per snippet
T_p = 30;   % number of samples preceding (and including) spike trough
thresh_factor = 4.5;    % xRMS threshold for detecting "spikes"

%% Create noise signal and look for spike-like snippets
noise = zeros(N,T_s);
count = 0;
while count < N
    %% Create noise signal
    % this is done multiple times because spike-like signals are rare and
    % making a single signal with enough "spikes" would be too large to
    % hold in memory. Instead we just keep making new noise signals until
    % we find enough "spikes".
    T = 10^6;   % length of noise signal (something long); MUST BE EVEN
    x = randn(1,T);
    X = fft(x);

    % prepare a vector for 1/f multiplication
    NumUniquePts = T/2 + 1;
    n = 1:NumUniquePts;
    n = sqrt(n);

    % multiplicate the left half of the spectrum so the power spectral density
    % is proportional to the frequency by factor 1/f, i.e. the
    % amplitudes are proportional to 1/sqrt(f)
    X(1:NumUniquePts) = X(1:NumUniquePts)./n.^2;

    % prepare a right half of the spectrum - a copy of the left one,
    % except the DC component and Nyquist frequency - they are unique
    X(NumUniquePts+1:T) = real(X(T/2:-1:2)) -1i*imag(X(T/2:-1:2));

    % prepare output vector y
    y = ifft(X);
    y = real(y(1, 1:T));

    % ensure unity standard deviation and zero mean value
    y = y - mean(y);
    yrms = sqrt(mean(y.^2));
    y = y/yrms;
    
    %% Filter noise signal like you would neural data
    hpFilt = designfilt('highpassiir', 'FilterOrder', 8, ...
                 'PassbandFrequency', 250, 'PassbandRipple', 0.2,...
                 'SampleRate', 10e3);
    y_filt = filtfilt(hpFilt,y);

    threshold = thresh_factor*rms(y_filt);   % threshold for spike detection
    crossed = false;
    i = T_s;
    while i<=T-T_s+T_p-1 && count<N
        %% Detect and extract snippets
        if ~crossed && y_filt(i)>threshold && y_filt(i-1)<threshold
            crossed = true; % signal has crossed threshold
        elseif crossed && y_filt(i)<y_filt(i-1)
            count = count+1;
            noise(count,:) = y_filt(i-T_p:i+T_s-T_p-1);
            crossed = false;
        end
        i = i+1;
    end
end
