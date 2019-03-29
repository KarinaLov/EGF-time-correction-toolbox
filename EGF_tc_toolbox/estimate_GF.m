function egfs = estimate_GF(settingsfile)
% Estimates the green's function from noise recorded on stationA and
% stationB between first_day and last_day
%
% Input:
%       settingsfile = text file where the input values are defined
%
% Output:
%       egfs = struct containing the estimated Green's function, lag time, 
%               the number of the days, and name of station pair
%
%
% Sub-function: read_settings.m, str2filename.m, gen_response.m, read_daily.m,
% prepros.m and cross_conv.m
%
% Written by Karina Løviknes 
% 

% Default values from settings file:
[network,stations,first_day,last_day,channels,location,num_stat_cc,Fq,filename,fileformat,pz_file,dateformat,deci,missingfiles,bpf,norm,wl,swl] = read_settings(settingsfile,'EGF');

validateattributes(stations,{'cell'},{'nonempty'});
nost = length(stations);

fd = datetime(first_day);
ld = datetime(last_day);
datevector = [fd:ld];
num_days = length(datevector); % Number of days
num_corr = num_days*24/swl;

tspd = Fq*60*60*24; % Total samples per day
lcc = 2*tspd/(24/wl)-1; % Length of cross correlation function

% Design the filter using the given cutoff frquencies and designfilt
df1 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',bpf(1),'HalfPowerFrequency2',bpf(2), ...
    'SampleRate',Fq,'DesignMethod','butter');

sp = 0; % Count the station pairs
ii = 0;
for jj = 1:nost-1
    
    % Loop over all the station pairs
    stationA = char(stations(jj))  
    
    % Pole zero file for station A: 
    pz_file_A = str2filename(pz_file,stationA,dateformat,'channels',channels,'network',network);
    respA = gen_response(tspd,Fq,pz_file_A);

    % Extract the data from the daily sac file (the sac-file needs to have
    % the format 'stationname-ch.yyyy-mm-dd.sac'
    SAdata = read_daily(network,stationA,channels,location,datevector,filename,fileformat,dateformat,Fq,deci,missingfiles);

for kk = 1:num_stat_cc-ii
     sp = sp+1;
    
    stationB = char(stations(jj+kk))
    
    pz_file_B=str2filename(pz_file,stationB,dateformat,'channels',channels,'network',network);
    respB = gen_response(tspd,Fq,pz_file_B);

    % Extract the data from all the daily files of station B
    SBdata = read_daily(network,stationB,channels,location,datevector,filename,fileformat,dateformat,Fq,deci,missingfiles);

    % Preallocate for speed:
    EGF = zeros(num_corr,lcc);
    
    pair = [stationA '-' stationB]
    dates = [char(first_day) '-' char(last_day)];
    
% ESTIMATE THE GREEN'S FUNCTION:
nk=24/swl;
for d = 1:num_days
    
    % Preprocess the data for each station, spesifying the
    % normalization steps to be applied:
    SAprosd = prepros(SAdata(d,:),Fq,df1,respA,norm);
    SBprosd = prepros(SBdata(d,:),Fq,df1,respB,norm);

    % Cross correlations:
    [EGF1 lag] = cross_conv(SAprosd,SBprosd,Fq,wl,swl);
    k=d*nk;
        
    EGF(k-(nk-1):k,:)=EGF1;
    
    % Save the daily cross correlations as SAC-files:
    %Header=struct('DELTA',1/Fq,'B',lag(1)/Fq,'E',lag(end)/Fq,'KSTNM',pair,'KHOLE',00,'KCMPNM',channels,'KNETWK',network,'NZDTTM',datevec(datevector(1)));
    %mksac(['Egf_' pair '_' datevector(d) '.SAC'],stack,datenum(first_day),Header)
    
end
stack=sum(EGF);

estimatedGF = struct('EGF',EGF,'lag',lag,'number_of_days',num_days,'pair',pair);
egfs(sp) = estimatedGF;
save(['Egf_' pair '_' dates '.mat'],'estimatedGF','-v7.3')

Header=struct('DELTA',1/Fq,'B',lag(1)/Fq,'E',lag(end)/Fq,'KSTNM',pair,'KHOLE',00,'KCMPNM',channels,'KNETWK',network,'NZDTTM',datevec(datevector(1)));
mksac(['Egf_' pair '_' dates '.SAC'],stack,datenum(first_day),Header)

end
% Make sure the stations are cross correlted with the rigth number of stations: 
if ii >= num_stat_cc
    ii = 0; 
else
    ii = ii + 1;
end
end
end




