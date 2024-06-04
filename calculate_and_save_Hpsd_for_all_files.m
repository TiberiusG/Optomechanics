% calculate and plot the power spectral density (PSD)
% of all .trc files containing the sequence written on line 16.
% Author: Tiberius Georgescu
% For this program to work, the functions 
% hpsd.m and ReadLeCroyBinaryWaveform should be in the same folder as 
% this matlab program

clear all;
close all;
clc;

f_min = 0;      % PSD display range min. frequency fmin
f_max = 2e5;    % PSD range fmax  
nos = 10;       % number of segments to which the time trace is divided and then averaged in the PSD estimate algorithm

files = dir('C4*.trc'); % search for all files of type .trc, containing the sequence 'C2'
nof = length(files);    % number of .trc files contatining 'C2'
%%
for i=1:nof
    file_name = files(i).name;  % get file name
%file_name = 'C1--Trace--00030.trc';

signal = ReadLeCroyBinaryWaveform(file_name);
% use the ReadLeCroyBinaryWaveform function to retrieve the recorded 
% oscilloscope information
N = length(signal.x); % extract the length of signal to use as input for psd function
time = signal.x;    % scope signal time coordinate
A = signal.y;       % scope signal amplitude or voltage coordinate
fs = signal.desc.fs;% sampling frequency used during recording

%% Plot the single sided PSD using the hpsd function, which must be in the same folder as this matlab program
[psd,f] = hpsd(A,nos,fs);

%%
% Reduce the PSD to the region where the peaks are
df= f(2)-f(1);                          % the frequency spacing
f1=f(f>= f_min -df/2 & f<= f_min+df/2 );% the frequency from the list that
                                        % is closest to our sought after value
index_f1=find(f==f1);                   % the index of the value of f1

                                        % identify values of f that are <= f_max
f2=f(f>= f_max -df/2 & f<=f_max+df/2);  % given the discretization df
index_f2=find(f==f2);                   % find the index of f2
f_reduced=f(index_f1:index_f2);         % create the reduced BW signal
psd_reduced=psd(index_f1:index_f2);     
% write the reduced psd to a .txt file
%dlmwrite([file_name(1:2),'_',file_name(end-5:end-4), '_f_PSD_',num2str(f_min),'_',num2str(f_max),'.txt'], [f_reduced psd_reduced], 'precision',15)
%dlmwrite([file_name(1:2),'_',file_name(end-5:end-4), '_f_PSD_',num2str(f_min),'_',num2str(f_max),'.txt'], [f_reduced psd_reduced], 'precision',15)

%% plot useful psd, around the peak of interest 
figure(2)
plot(f_reduced/1e3, psd_reduced)        % plot single sided PSD versus frequency
set(gca, 'Yscale', 'log')               % set y scale to log10, to better see the peaks
axis([f_min/1e3 f_max/1e3 -inf inf])    % confine the displayed values to these intervals


title([file_name(1:2),' ', file_name(end-5:end-4),' PSD between ',num2str(f_min/1e3),'-',num2str(f_max/1e3),'kHz'])
xlabel('Frequency [kHz]')
ylabel('PSD [V^2/Hz] ')
saveas(gcf, [get(get(gca,'title'),'string'),'.png'])    % save as .fig
saveas(gcf, get(get(gca,'title'),'string'))             % save as .png
close(gcf) % close figure after plotting
end

