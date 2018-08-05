%% Script for Detecting SPW-R events that last for at least 15ms

%% Directory
cd 'x:\03. Lab Procedures and Protocols\MATLABToolbox\chronux\spectral_analysis\continuous';

%% LFP Load
% load '';
% 

%% Build the Gaussian Kernel (4-ms s.d - per Karlsson, 2009)
pts_ms=params.Fs/1000; %points/ms in ms scale
pts_env=(1:(ceil(pts_ms*15))); pts_env=ones(1,length(pts_env)); %num data points expected within a 15ms window
sig=pts_ms*4; %4ms worth of data (per Karlsson, 2009 s.d.)
max=ceil(3*sig);
gaus=fspecial('gaussian',2*max+1,sig); %generate 2D gaussian
n=ceil(length(gaus)/2);%round to the half-way row
gaus=gaus(n,:); %make vector equal to middle-row
gaus=gaus./sum(gaus); %correct gaussian kernel to sum to 1
%% Pull out LFP from given channel
data=CSC(n).whole; %pick the channel

%% Bandpass filter
data_f=skaggs_filter_var(data,140,220,params.Fs); %filters for 150-250Hz

%% Hilbert Transform
data_h=abs(hilbert(data_f));

%% Convolve with Gaussian Kernel - smoothing
data_s=conv(data_h,gaus,'same');

%% standardize amplitude
data_s_z=zscore(data_s);

%% Identify + 3stdv from mean abs. amplitude
data_amp_avg=mean(data_s_z); %calculate the avg amplitude
data_std=std(data_s_z); %standard dev for absolute amplitudes
data_normax=data_amp_avg+(data_std*3); %maximum expected value in normal distribution defined as 3stds above mean


%% SPW-R detection
data_flag=data_s_z>data_normax;% returns boolean for all cells in the absoluted signal that are greater than normal max
data_env=strfind(data_flag,pts_env);% cases in which minimum length of TRUE is found in boolean signal
[data_flagCol, data_reps, data_ind]=RunLength(data_flag);% Collapsed unbroken sequence of same TRUE/FALSE value; number of items in collapse; cell column value in original signal in which TRUE/FALSE value change
data_multiples = find(data_reps>length(pts_env));%cases in which the number of collapsed values is greater than minimum
data_SWRchk=data_ind(data_multiples); %the first of cell in which the collapsed sequence started

for q=1:length(data_SWRchk);
    if data_flag(data_SWRchk(q))==1; %if the sequence met the duration criteria, is that a sequence of 1's
        data_SWRstart(q)=data_SWRchk(q); %if so, enter here
    else
        data_SWRstart(q)=NaN; %if not enter as NaN
    end
end

data_SWRstart(isnan(data_SWRstart)) = []; %clear NaN's

data_SWRcount=length(data_SWRstart); %total number, for each trial, of SPW-R events

%% Visualize SWRs
%extract original signal for each event
for d = 1:data_SWRcount; %for the number of detected events
    SWRs(d,:)=data(1,data_SWRstart(1,d)-250:data_SWRstart(1,d)+250); %make a matrix with rows representing the event, and columns populated with signal values from start-250 to start+250
end

%filter for ripple band
for d=1:data_SWRcount;
    SWRs_f(d,:)=skaggs_filter_var(SWRs(d,1:end),140,220,params.Fs);
end

%filter for theta band
for d=1:data_SWRcount;
    SWRs_f_the(d,:)=skaggs_filter_var(SWRs(d,1:end),6,12,params.Fs);
end

%plot it
subplot 311; plot(SWRs(513,1:end)); subplot 312; plot(SWRs_f_the(513,1:end)); subplot 313; plot(SWRs_f(513,1:end)); %using event #200 as an example

%% Save output
save ('X:\08. Lab personnel\Current\David\Projects\Re Suppression - HC Modulation\2. Output\Ephys\Data Analysis\Santiago\DelayPhaseLtNLt_SWRs.mat','-v7.3');
