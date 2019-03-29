function [fdelay fdelayc] = invert_TD(settingsfile,varargin)
% The delays found for each station pair is inverted to find the final
% time delay of each station
%
% Input:
%       settingsfile = text file where the input values are defined
%
% Output:
%     fdelay = the final delay after inversion
%     fdelayc = the final delay after inversion and correction
%
%
% Sub-function: read_settings.m
%
% Written by Karina Løviknes
%

% Default values from settings file
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,swl,RCS,yaxis] = read_settings(settingsfile,'INVERT');

validateattributes(stations,{'cell'},{'nonempty'});
nost=length(stations);

fd = datetime(first_day);
ld = datetime(last_day);
datevector = [fd:ld];
num_days = length(datevector); % Number of days

nsp = nost*num_stat_cc/2; % Number of station pair 

sp = 0; % Count the stationpairs
G = zeros(nsp,nost); % Matrix for inversion:

for jj = 1:nost
    % Loop over all the station pairs
    stationA = char(stations(jj));    

for kk = 1:num_stat_cc-(jj-1)
    sp = sp+1;
    
    % Fill in the inversion matrix with the station pair combinations:
    G(sp,jj) = 1; 
    G(sp,kk+jj) = -1;
    
    stationB = char(stations(jj+kk));

    pair = [stationA '-' stationB];
    dates = [char(first_day) '-' char(last_day)];
    
    % Extract measured time delays from file
    filename = ['TD_' pair '_' dates '.mat']; 
    if exist(filename,'file') % Check that the file exists
        file = load(filename);
        timedelay(sp,:) = file.timedelay.timedelay;
        timedelay0(sp,:) = file.timedelay.timedelay0;
        num_days = file.timedelay.number_of_days;
        num_corr = length(timedelay(sp,:));       
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: TD_' pair '_' dates '.mat' ])
    end
end
end

% Inversion matrix:
G

% Calculate and plot the final time delay using the pseudo inverse
fdelay = pinv(G)*timedelay;

% Correct the relative time delay by setting the most realiable station to zero:
fdelay1 = fdelay;
if isempty(RCS)
    % Find the station with the smallest number of measured time delays
    for i = 1:nost
        fdelay1(i,find(isnan(fdelay1(i,:)))) = 0;
        sum_m(i,:) = sum(abs(fdelay1(i,:)));
    end

    realm = find(sum_m==min(sum_m));
    
else
    realm = find(strcmp(stations,cellstr(RCS)));
end

% Correct to absolute time delay by removing the deay of the reliable
% station:
for i = 1:nost
    fdelayc(i,:) = fdelay(i,:)-fdelay(realm,:);
end
    
if strcmp(varargin,'plot')
    % Plot the corrected inverted timedelays 
    figure
    for ff = 1:nost
        subplot(nost,1,ff)
        plot(fdelayc(ff,:)/Fq)
        axis([0 num_corr yaxis])
        title(['Corrected timedelay for ' char(stations(ff))])
        xlabel('days')
        hold on
    end
    hold off 
end
end