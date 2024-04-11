load MATLAB\DATA\Gait\30HZ\30hz_ST.mat;


% constants definition
emg_data1=emg_data(40000:100000);
samplingRate = 2500;  %see Biometrics documentation
Fnyq = samplingRate/3;
length_sec = size(emg_data1,1)/samplingRate;
Ts = 1/samplingRate;
ii = Ts:Ts:length_sec;


% zero-phase filtering
Rectified_data = filtfilt(ones(1,500)/500,1,abs(emg_data1));


% Low pass filtering
N = 4;              % filter order
f_cutOff = 24;      % cut-off frequency for low-pass filter
[D1, C1]  = butter(N, f_cutOff*1.116/Fnyq,'low'); % Butterworth filter
Filtered_data = filtfilt(D1, C1, Rectified_data);


figure;
hold on
plot(ii, emg_data1, 'blue');
%plot(t,factor(i)*Filtered_data(:,i), 'red'); -> if you want to plot the
%filtered data from the square of the raw data
plot(ii,Filtered_data, 'red');
xlabel('Time (s)');
ylabel('Signal (mV)');
legend('raw signal','Filtered signal (from absolute value of raw data)');


figure;
subplot(2,1,1);
plot(Rectified_data); hold on;
subplot(2,1,2);
plot(Filtered_data);  hold on;

% threshold = 0.001;
% ii=1;
%  while abs(Filtered_data(ii+20,1)-Filtered_data(ii,1))>0.001     %avoid to take into account the little slope at the beginning of the measurement
%        ii = ii+1;
%  end
%         while (Filtered_data(ii+50,1)-Filtered_data(ii,1))/samplingRate<threshold
%         ii= ii+1;
%         end
%         init_ext = ii;
% 
% marg = 200;
% 
%    Filtered_data = Filtered_data(init_ext-marg:end);
%    Time_ext = size(Filtered_data,1); 
%     
% % Time normalization
% xq = [1:min(Time_ext)];
% if size(Filtered_data,1)~=xq
%     Ext_timeNorm = interp1(Filtered_data,xq);
% else
%     Ext_timeNorm = Filtered_data;
% end
% Time_ext_pc = 1/Time_ext*100:1/Time_ext*100:100;
% 
% 
% % for extension
% Ext_mean = [];
% Ext_std = [];
% EMG = Ext_timeNorm;
%     Ext_mean=[Ext_mean mean(EMG,2)];
%     Ext_std=[Ext_std std(EMG,0,2)];
% 
% 
%     Ext_mean_normalized = Ext_mean/(max(Ext_mean));
% 
% 
%     figure();
% plot(Time_ext_pc,Ext_mean_normalized);
% xlabel('cycle(%)');
% ylabel('EMG normalized(% of maximum value)');
% title('EMG normalized for knee extension');
% 
% 
