% open trc file, calculate the power spectral density (PSD), plot roi, fit to a peak
% Author: Tiberius Georgescu

clear all;          % close all matlab files and delete all variable values
close all;
clc;

%Variables
tic;            % start counting time
f_min = 102e3;%71e3;%36e3;   % min. freq. of the fit interval, CHANGE MANUALLY FOR EACH PEAK 
f_max = 112e3;%81e3;%57e3;   % max. freq., CHANGE MANUALLY FOR EACH PEAK  
f0 = 107e3;%76e3;%46.5e3;    % peak frequency, CHANGE MANUALLY FOR EACH PEAK 
bw = f_max - f_min; % bandwidth
nos = 10;       % number of segments to which the time trace is divided and then averaged in the PSD estimate algorithm
k = 3;          % the peak name index 123 = zyx, CHANGE MANUALLY FOR EACH PEAK 

%open .trc file and extract signal
file_name = 'C4--XX--00007.trc';  % the name of the oscilloscope .trc file, ADD MANUALLY 
%            1234567890123456789012 
signal=ReadLeCroyBinaryWaveform(file_name); % get the data from the .trc file, using the function 'ReadLeCroyBinaryWaveform' 
%file_name = scope_file_name;

N = length(signal.x);   % extract the length of signal to use as input for psd function
time = signal.x;        % scope time record
A = signal.y;           % scope voltage signal
fs = signal.desc.fs;    % scope sampling frequency used for acquisition

%Calculate the PSD 
[psd,f] = hpsd(A,nos,fs);
% Reduce the PSD to the region where the peaks are
df= f(2)-f(1);                          % the frequency spacing
f1=f(f>= f_min -df/2 & f<= f_min+df/2 );% the frequency from the list that
                                        % is closest to f_min
index_f1=find(f==f1);                   % the index of the value of f1

                                        % identify values of f that are <= f_max
f2=f(f>= f_max -df/2 & f<=f_max+df/2);  % given the discretization df
index_f2=find(f==f2);                   % find the index of f2
f_useful=f(index_f1:index_f2);          % create the reduced BW signal
psd_useful=psd(index_f1:index_f2);      % the reduced BW PSD
% write the reduced psd to a .txt file
% dlmwrite([file_name(1:8), '_f_PSD_',num2str(f_min),'_',num2str(f_max),'.txt'], [f_reduced psd_reduced], 'precision',15)

% plot useful psd, around the peak of interest 
% figure(2)
% plot(f_useful/1e3, psd_useful)
% set(gca, 'Yscale', 'log')
% title([file_name(1:2),' ', file_name(end-5:end-4),' PSD between ',num2str(f_min/1e3),'-',num2str(f_max/1e3),'kHz'])
% xlabel('Frequency [kHz]')
% ylabel('PSD [V^2/Hz] ')

% find the maximum
Ap = max(psd_useful);       % Amplitude of peak = max value of reduced PSD
loc = find(psd_useful==Ap); % find index of max value
fp = f_useful(loc);         % get the frequency associated to that index

% start subplot loop
pk_names = '123456789';     % peak names
pkname = pk_names(k);       % get the k character from pk_names string

figure(3)                   % plot the reduced/useful PSD
plot(f_useful/1e3, psd_useful)% divide freq. by 1e3 to plot in kHz
set(gca, 'Yscale','log')    % set scale of y axis to log10 
title(['PSD z peak ',num2str(bw/1e3),' kHz BW']) % include in title PSD BW
xlabel('Frequency [kHz]')   % label of x axis
ylabel(' V^2/Hz ')          % label of y axis

%% ESTIMATE FITTING PARAMETERS FROM THE EXPERIMENTAL DATA SET, AS YOU WOULD DO BY EYE

hm_t = Ap / 2; % calculate half of the peak height

dif = zeros(length(psd_useful),1); %create an array that will store the differences
%* Calculate the difference between each Y point and Ymax/2, with sign
for i = 1 : length(psd_useful)
    dif(i) = psd_useful(i) - hm_t;
end
%* Find the zero crossings and select from the two points, dy<0 & dy>0,
%  the one that is closest to Ymax/2 
points = zeros(length(psd_useful),1);%points will become 1 wherever y = Ymax/2
for j = 1 : length(psd_useful) - 1  %if the last point of the psd, which is not included in this count, is zero that does not affect the result
    if dif(j) == 0                          %if Yj - Ymax/2 = 0,
        points(j) = 1;                      %then points(j) = 1, which means that at position j is one of the points I am after
    elseif dif(j) < 0 && dif(j+1) > 0 ...   %if the Ymax/2 is between j and j+1
        && abs(dif(j)) <= abs(dif(j+1))     %and j is closer to Ymax/2 or j and j+1 are at exactly the same distance from max/2,
    points(j) = 1;                          %then pick the j position as the correct point
                                            
    elseif  dif(j) < 0 && dif(j+1) > 0 ...  
        && abs(dif(j)) > abs(dif(j+1))
        points(j+1)=1;
    
    elseif dif(j) > 0 && dif(j+1) < 0 ...
        && abs(dif(j)) <= abs(dif(j+1))
        points(j) = 1;
    elseif dif(j) > 0 && dif(j+1) < 0 ...
        && abs(dif(j)) > abs(dif(j+1))
        points(j+1) = 1;
    end
end
%* find the number of points different from zero in points
np = length(points(points==1)); % nr of elements equal to 1 from 'points'
%* Now find the frequencies of the y = Ymax/2 points
x_p = points.*f_useful;     %the x and y coords of all the points whose 
y_p = points.*psd_useful;   %y(point) = Ymax/2
%fwhm Gamma = fmax - f min, so we are just interested in the x comp.
i1 = min(find(x_p~=0)); % index of lowest freq. whose height is Ap/2
i2 = max(find(x_p~=0)); % index of highest freq. whose height is Ap/2

%% The initial fit coefficients
g = f_useful(i2) - f_useful(i1);    % damping rate, or width at half maximum
A1 = sum(psd_useful)*df;            % the area under the PSD, found by integration
n = (psd_useful(1)+psd_useful(end));% the noise floor

x = f_useful;      % x coordinate
y = psd_useful;    % y coord.
 
S = A1/pi * 2*f0^2*g ./...              % the fit function constructed 
    ( (x.^2 - f0^2).^2 + x.^2*g^2 ) + n;% with the initial coefficients
 

figure(4)
hold on;
scatter(x/1e3, y)
set(gca, 'Yscale','log')
title([file_name(1:2),' ', file_name(end-5:end-4),' PSD fit ',pkname])
xlabel('frequency [kHz]')
ylabel(' PSD [V^2s/rad]')
%figure(8)
plot(x/1e3,S,'color','red','LineWidth',1)%,'--')
legend(' Exp','Fit')
hold off;
saveas(gcf, [get(get(gca,'title'),'string'),'.png'])

%% Write the preliminary fit coefficients with 15 decimal points precision
% out_ref_1 = [A1
%     f0
%     g ];
% S_data_out = [f_useful psd_useful S ];
% dlmwrite([file_name(1:2),' ', file_name(end-5:end-4), '_h_f_PSD_fit_',pkname,'.txt'], S_data_out, 'precision',15);
% dlmwrite([file_name(1:2),' ', file_name(end-5:end-4), '_h_A1f0gn_out_ref_',pkname,'.txt'], out_ref_1, 'precision', 15);

%% * the fitting algorithm
tol = 1e-6; % the tolerance

err=inf;    % initiate error as infinite

%[y]=[f]+[Z][da]+[e]
while err>tol
%compute column vector f, the actual fitting function
fit = A1/pi * 2*f0^2*g ./...
    ( (x.^2 - f0^2).^2 + x.^2*g^2 ) + n;
%compute the column da =w measured- fitted
d=y-fit;
%construct the Z matrix, with partial derivatives
%   z1 = dfit/dA1
z1 = (2*f0^2*g)./(pi*((x.^2 - f0^2).^2 + x.^2*g.^2));
%   z2 = dfit/df0
z2 = (4*A1*f0*g)./(pi*((x.^2 - f0^2).^2 + x.^2*g^2)) + (8*A1*f0^3*g*(x.^2 - f0^2))./(pi*((x.^2 - f0^2).^2 + x.^2*g^2).^2);
% z3 = dfit/dg
z3 = (2*A1*f0^2)./(pi*((x.^2 - f0^2).^2 + x.^2*g^2)) - (4*A1*x.^2*f0^2*g^2)./(pi*((x.^2 - f0^2).^2 + x.^2*g^2).^2);
% z4 = dfit/dn 
z4 = ones(length(fit),1);
Z=[z1 z2 z3 z4];
%UPDATE COEFFICIENTS
%da is the error
da = ( Z.' * Z ) \ ( Z.' * d );

%a = a + da(1);
A1 = A1 + da(1);
%b = b + da(2);
f0 = f0 + da(2);
%c = c + da(3);
g = g+ da(3);
n = n + da(4);

relative_err= da./[A1; f0; g; n];
err=max( abs( da./[A1; f0; g; n]));
% By including the figure in the While loop, we can see how the fit evolved     
% figure(5)
% hold on;
% plot(x/1e3, y)
% plot(x/1e3,fit,'color','red','LineWidth',1)
% set(gca,'Yscale','log')
% title([file_name(14:16),' Fitted function VS experimental data PSD'])
% xlabel('Frequency [kHz]')
% ylabel('Amplitude [V^2/Hz]')
% legend([file_name(14:16),' Experimental data'],'Fitted function')
% saveas(gcf, get(get(gca,'title'),'string'))
% %saveas(gcf, [file_name(7:9),get(get(gca,'title'),'string')])
% hold off;
%err=0;
end
%% After fit figure
figure(6)
hold on;    % hold on to plot both PSD and fit on the same graph
scatter(x/1e3, y)
plot(x/1e3,fit,'color','red','LineWidth',1)
set(gca,'Yscale','log')
title([file_name(1:2),' ', file_name(end-5:end-4),' Fitted function VS experimental data PSD ',pkname])
xlabel('Frequency [kHz]')
ylabel('Amplitude [V^2/Hz]')
legend([file_name(1:2),' ', file_name(end-5:end-4),' Experimental data'],'Fitted function')
saveas(gcf, [get(get(gca,'title'),'string'),'.png'])    % save figure as .png, with the name given by the figure title
%saveas(gcf, [file_name(7:9),get(get(gca,'title'),'string')])% save figure
%as .fig 
hold off;
%% Make an output array containing the values of the found fit coefficients
out_ref = [A1 
    f0
    g
    n
    ];

data_out = [f_useful  psd_useful fit];  % 2 col. matrix containing freq. and PSD
% save fit coefs. in .txt file
dlmwrite([file_name(1:2),' ', file_name(end-5:end-4), '_A1f0gn_out_ref_',pkname,'.txt'], out_ref, 'precision', 15);
% save freq. and psd in a .txt file
dlmwrite([file_name(1:2),' ', file_name(end-5:end-4), '_f_PSD_fit_',pkname,'.txt'], data_out, 'precision', 15);
% the 15 digit precision is essential when saving in .txt file; the default
% precision is too small and will cause problems when looking at data

toc; % stop timer.
% the difference between tic and toc is the time it took the computer to
% perform all computations

% references:
% fit algorithm:
% https://www.youtube.com/watch?v=7V9KP7BkO1M&ab_channel=EMPossible
% 2017_Hebestreit_Thesis_Thermal_properties_of_levitated_particles
