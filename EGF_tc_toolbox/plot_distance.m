function h = plot_distance(settingsfile) 
% Calculates the inter-station pair distance and plots the stacked noise 
% correlations with distance
%
% Input: 
%       settingsfile = textfile containin where the input values are defined 
%
% Output:
%       h = figuren
%
%
% Sub-function: read_settings.m, str2filename.m, read_info_pz.m and
%       calc_dist.m
%
% Written by Karina Løviknes 
% 

% Read all the default values from the settings file
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,xaxis,pz_file,dateformat,lag_red] = read_settings(settingsfile,'DISTANCE');
    
validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);

% Loop over all the station pairs:
sp=0; % Count the station pairs
h = figure;
for jj=1:nost-1
    stationA=char(stations(jj)); 
    
    % Read coordinates for station A from the pole zero file 
    pz_file_A = str2filename(pz_file,stationA,dateformat,'channels',channels,'network',network);
    [longA,latA] = read_info_pz(pz_file_A,'Coordinates');

for kk=1:num_stat_cc-(jj-1)
     sp=sp+1;
    
    stationB=char(stations(jj+kk));

    pair=[stationA '-' stationB];
    dates=[char(first_day) '-' char(last_day)];
    
    filename1=['Egf_' pair '_' dates '.mat'];
    if exist(filename1,'file') % Check that the file exists
        file1=load(filename1);
        EGF=file1.estimatedGF.EGF;
        lag=file1.estimatedGF.lag;
        num_days = file1.estimatedGF.number_of_days;
        num_corr = length(EGF(:,1));
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: Egf_' pair '_' dates '.mat' ])
    end

    filename2 = ['TD_' pair '_' dates '.mat']; 
    if exist(filename2,'file') % Check that the file exists
        file2 = load(filename2);
        timedelay = file2.timedelay.timedelay;
        timedelay0 = file2.timedelay.timedelay0;
        ref(sp,:) = file2.timedelay.reference;      
    else
        error(['Cannot find a mat.file with an estimated greens function for stationpair ' pair '. Fileformat must be: TD_' pair '_' dates '.mat' ])
    end
          
    % Make lag-vector
    lr = length(ref(sp,:))-1;
    lag = -lr/2:lr/2;
    
    % Read coordinates for station B from the pole zero file 
    pz_file_B = str2filename(pz_file,stationB,dateformat,'channels',channels,'network',network);
    [longB,latB] = read_info_pz(pz_file_B,'Coordinates');
    
    % Calculate distance:
    in_st_dist(sp) = calc_dist(latA,longA,latB,longB);

    refd = (ref(sp,:)/max(abs(ref(sp,:))))*10+in_st_dist(sp);
    plot(lag/Fq,refd)
    axis([xaxis 0 max(in_st_dist)+10])
    title('Cross correlations in stationpair distance')
    xlabel('Time lag (s)'), ylabel('Distance (km)')
    hold on  
end
end
hold off
end



