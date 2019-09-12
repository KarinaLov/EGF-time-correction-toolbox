function h = plot_distance(settingsfile) 
% Calculates the inter-station pair distance and plots the stacked noise 
% correlations with distance
%
% Input: 
%       settingsfile = textfile containin where the input values are defined 
%
% Output:
%       h = figure
%
%
% Sub-function: read_settings.m, str2filename.m, read_info_pz.m and
%       calc_dist.m
%
% Written by Karina Løviknes 
% 

% Read all the default values from the settings file
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,xaxis,pz_file,dateformat,lag_red,bpf,dates] = read_settings(settingsfile,'DISTANCE');
    
validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);
nch = length(channels);
nsp = nch*nost*num_stat_cc/2; % Number of station pair 

dates1 = [char(first_day) '-' char(last_day)];
if isempty(dates)
    dates2 = dates1;
else
    dates2 = [char(dates(1)) '-' char(dates(2))];
end

% If spesified, filter the daily cross correlations:
filttxt = '';
if ~isempty(bpf)
    % Design the filter using the given cutoff frquencies and designfilt
    ord = 4; % Filter order, default
    df = designfilt('bandpassiir','FilterOrder',ord, ...
        'HalfPowerFrequency1',bpf(1),'HalfPowerFrequency2',bpf(2), ...
        'SampleRate',Fq,'DesignMethod','butter');
    filttxt = [' - filt: ' num2str(bpf(1)) '-' num2str(bpf(2))];
end

% Preallocate for speed
in_st_dist=zeros(1,nsp);

% Loop over all the station pairs:
sp = 0; % Count the station pairs
ii = 0;
for ch = 1:nch
    channel = channels(ch);
for jj=1:nost
    stationA=char(stations(jj)); 
    
    % Read coordinates for station A from the pole zero file 
    pz_file_A = str2filename(pz_file,stationA,dateformat,'channels',channels,'network',network);
    [longA,latA] = read_info_pz(pz_file_A,'Coordinates');

for kk=1:num_stat_cc-ii
     sp=sp+1;
    
    stationB = char(stations(jj+kk));

    pair = [stationA '-' stationB '-' channel]
    %pair = [stationA '-' stationB]
    pairs(sp) = {pair};

    filename2 = ['TD_' pair '_' dates2 '.mat']; 
    if exist(filename2,'file') % Check that the file exists
        file2 = load(filename2);
        timedelay = file2.timedelay.timedelay;
        timedelay0 = file2.timedelay.timedelay0;
        ref = file2.timedelay.reference;      
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: TD_' pair '_' dates2 '.mat' ])
    end
          
    % Make lag-vector
    lr = length(ref)-1;
    lag = -lr/2:lr/2;
    
    % If spesified, filter the reference
    if ~isempty(bpf)
        % Filter the refrence
        reff(sp,:) = filtfilt(df,ref);
    else
        reff(sp,:) = ref;
    end
    
    % Read coordinates for station B from the pole zero file 
    pz_file_B = str2filename(pz_file,stationB,dateformat,'channels',channel,'network',network);
    [longB,latB] = read_info_pz(pz_file_B,'Coordinates');
    
    % Calculate distance:
    in_st_dist(sp) = calc_dist(latA,longA,latB,longB)*1000;
end
% Make sure the stations are cross correlted with the rigth number of stations: 
if ii >= num_stat_cc
    ii = 0; 
else
    ii = ii + 1;
end
end
end

maxdist = max(in_st_dist)
btwdist = maxdist/sp;
mindist = min(in_st_dist);
in_st_dist

h = figure;
for ii=1:sp
    refd = (reff(ii,:)/max(abs(reff(ii,:))))*(10)+in_st_dist(ii);
    plot(lag/Fq,refd)
    axis([xaxis mindist-15 maxdist+15])
    title(['Cross correlations in stationpair distance' filttxt])
    xlabel('Time lag (s)'), ylabel('Distance (m)')
    text(3.1,in_st_dist(ii)+(mindist/sp),pairs(ii))
    hold on  
end
hold off
end




