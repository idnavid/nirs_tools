function [fnirs_signal,hrf_true,X] = simulate_nirs_data(length_signal,length_hrf,...
                                                        pulse_on_width,pulse_off_width,fs)
% Create Haemodynamic response to pulse train stimulus. 
% Inputs:
%       length_signal   : length of output fNIRS signal (samples)
%       length_hrf      : length of haemodynamic response function (HRF) (samples)
%       pulse_on_width  : length during which the stimulus is on (high)
%       pulse_off_width : length during which the stimulus is off (low) (samples)
%       fs              : sampling frequency
%       
% Outputs:
%       fnirs_signal    : simulated fnris signal
%       hrf_true        : HRF
%       X               : autoregressive basis vectors (fnirs_signal = X*hrf_true + unwanted_signal)

hrf_true = hrf(length_signal,length_hrf);

t = (0:1:length_signal- 1)/fs;   %Time vector 
delay = (0:(pulse_on_width+pulse_off_width):length_signal-1)/fs; %delay vector
pulse_train = pulstran(t,delay,'rectpuls',pulse_on_width/fs); 

R = [pulse_train(1) zeros(1,length_hrf-1)];
X = toeplitz(pulse_train,R);
fnirs_signal = X*hrf_true';
fnirs_signal = fnirs_signal(1:length_signal);
drift_signal = 0.5*drift_mri(length_signal);
fnirs_signal = fnirs_signal(:)+drift_signal(:);
