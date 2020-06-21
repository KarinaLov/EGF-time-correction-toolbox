function [network, stations, first_day, last_day, channels, location,...
    num_stat_cc, Fq, varargout] = read_settings(settingsfile, varargin)
% Reads the settings file
%
% Input:
%       settingsfile=the name of the settingsfile
%       varargin: spesfification about the type of values to read
%               'EGF': values used to estimate the Green's funcion
%               'TD': values used to measure the time shift
%               'INVERT': values used to invert the time shift
%               'PLOT': values used to plot the estimated Green's function
%                       and/or time shifts
%               'FILTER': values used to filter and plot the Green's
%                       function
%               'DISTANCE': values used to plot the Green's function with
%                       distance
%               'CORRECT': values used to correct the data for the measured
%                       the time shift
%
% Output:
%       Fq=sampling frequency
%       filetype=spesification about type of input file, either sac or
%                   miniseed
%       num_stat_cc=Number of stations each station is cross correlated with
%       varargout
%
% Written by Karina Loviknes 
% 

% The optional variables are empty or default if not specified
vouts = {};
vouts{1} = 'network.stationname.00.HHchannels.D.yyyy.ddd.000000.SAC';
vouts{2} = 'sac';
vouts{3} = 'SAC_PZs_network_stationname_HHchannels.pz';
vouts{4} = 'yyyy-mm-dd';
vouts{6} = 6; % missingfiles tolerance
vouts{9} = 24; % cross correlation window length
vouts{10} = 24; % cross correlation stacking length
vouts{11} = 100; % overlapping percent cross correlation
datesm = {};
bpfm = {};
iterations = 3;
lag_red = 5000;
stackperiod = split('whole 0');
thr = 0.4;
fit = 'direct';
fitperiod = {};
RCS = {};
xaxis = [-150 150];
yaxis = [-100 100];
titl = 'name';
bpfp = {};
fc1 = {};
fc2 = {};

% Read the textfile into a table:
txt = readtable(settingsfile,'Delimiter','=', 'CommentStyle',{'%'},'ReadVariableNames',false);

szt=size(txt);
nc=szt(1); % Number of columms

vo=1; % Count the number of varargout variables
 for j=1:nc
    var_name=char(txt{j,1})
    var=char(txt{j,2})
    
    % Extract the variables:
    if strcmp(var_name, 'network')
        network=var
    
    elseif strcmp(var_name, 'stations')
        stations=split(var);
    
    elseif strcmp(var_name, 'first_day')
        first_day=var;
        
    elseif strcmp(var_name, 'last_day')
        last_day=var;
        
    elseif strcmp(var_name, 'channels')
        channels=var;
         
    elseif strcmp(var_name, 'location')
        location=var;

    elseif strcmp(var_name, 'num_stat_cc')
        % Number of station each station is cross correlated with
        num_stat_cc=str2num(var);  
            
    elseif strcmp(var_name, 'Fq')
        % Sampling rate
        Fq=str2num(var);
        
    elseif strcmp(var_name, 'filename')
        % Filename format
        vouts{1}=var;
        
    elseif strcmp(var_name, 'fileformat')
        % Fileformat
        vouts{2}=var;
        
    elseif strcmp(var_name, 'pz_file')
        % Pole zero file filename format
        vouts{3}=var;
        
    elseif strcmp(var_name, 'dateformat')
        % Format of the date in the filename
        vouts{4}=var;
        
    elseif strcmp(var_name, 'decimate')
        % How much to downsample 
        vouts{5}=str2num(var);  
            
    elseif strcmp(var_name, 'missingfiles')
        % Number of missing files in a row tolerated
        vouts{6}=str2num(var);  
        
    elseif strcmp(var_name, 'bpf')
        % Bandpass filter to apply during pre processing 
        vouts{7}=str2num(var); 

    elseif strcmp(var_name, 'norm')
        % Normalization during pre processing
        vouts{8}=split(var);
        
    elseif strcmp(var_name, 'wl')
        % Window length for cross correlation
        vouts{9}=str2num(var);
            
    elseif strcmp(var_name, 'swl')
        % Stacking widow length
        vouts{10}=str2num(var);
        
    elseif strcmp(var_name, 'perco')
        % Stacking widow length
        vouts{11}=str2num(var);

    elseif strcmp(var_name, 'datesm')
        % Bandpass filter to apply before measuring time errors
        datesm=split(var);
        
    elseif strcmp(var_name, 'bpfm')
        % Bandpass filter to apply before measuring time errors
        bpfm=str2num(var); 

    elseif strcmp(var_name, 'iterations')
        % Number of iterations to run the measuring process
        iterations=str2num(var); 
        
    elseif strcmp(var_name, 'lag_red')
        % Time lag to use for measuring time errors and plotting 
        lag_red=str2num(var); 
        
    elseif strcmp(var_name, 'stackperiod')
        % Time period to stack over for to use for reference trace during
        % timing error measuring 
        stackperiod=split(var);   
            
    elseif strcmp(var_name, 'signalpart')
        % Time period to stack over for to use for reference trace during
        % timing error measuring 
        signal_part = var; 
        
    elseif strcmp(var_name, 'threshold')
        % Time period to stack over for to use for reference trace during
        % timing error measuring 
        thr = str2num(var); 
        
    elseif strcmp(var_name, 'fit')
        % Station with reliable clock for the invertion 
        fit = var;  

    elseif strcmp(var_name, 'fitperiod')
        % Station with reliable clock for the invertion 
        fitperiod = split(var);  
            
    elseif strcmp(var_name, 'reference_clock_station')
        % Station with reliable clock for the invertion 
        RCS = var;   

    elseif strcmp(var_name, 'xaxis')
        % X-axis for the EGF
        xaxis=str2num(var); 

    elseif strcmp(var_name, 'yaxis')
        % Y-axsis for the measured time shift
        yaxis=str2num(var); 
        
    elseif strcmp(var_name, 'titl')
        % PLot title
        titl=var;
        
    elseif strcmp(var_name, 'bpfp')
        % Bandpass filter before plot
        bpfp=str2num(var);
        
    elseif strcmp(var_name, 'cutoff_freq1')
        % Bandpass filter before plot
        fc1 = str2num(var);
        
    elseif strcmp(var_name, 'cutoff_freq2')
        % Bandpass filter before plot
        fc2 = str2num(var);
            
    elseif strcmp(var_name, 'filenameO')
        % Filename format
        filenameO=var;
        
    elseif strcmp(var_name, 'fileformatO')
        % Fileformat
        fileformatO=var;
        
    elseif strcmp(var_name, 'dateformatO')
        % Format of the date in the filename
        dateformatO=var;
    end       
 end

vo = length(vouts); % number of possible varargout 
if isempty(varargin)
    % Read the entire file
    varargout=vouts;
    varargout{vo+1}=datesm;
    varargout{vo+2}=bpfm;
    varargout{vo+3}=iterations;
    varargout{vo+4}=stackperiod;
    varargout{vo+5}=fit;
    varargout{vo+6}=fitperiod;
    varargout{vo+7}=RCS;
    varargout{vo+8}=xaxis;
    varargout{vo+9}=yaxis;
    varargout{vo+10}=titl;
    varargout{vo+11}=bpfp;
    varargout{vo+12}=lag_red;
    varargout{vo+13}=fc1;
    varargout{vo+14}=fc2;
    
elseif strcmp(varargin{1},'EGF')
    % Read the variables used to estimate the Green's function
    varargout=vouts;

elseif strcmp(varargin{1},'TD')
    % Read the variables used to measure timeshifts
    varargout{1}=datesm;
    varargout{2}=bpfm;
    varargout{3}=iterations;
    varargout{4}=lag_red;
    varargout{5}=stackperiod;
    varargout{6}=signal_part;
    varargout{7}=thr;
    varargout{8}=fit;
    varargout{9}=fitperiod;
    
elseif strcmp(varargin{1},'INVERT')
    % Read the variables used to invert
    varargout{1}=fit;
    varargout{2}=fitperiod;
    varargout{3}=RCS;
    varargout{4}=yaxis;
    varargout{5}=datesm;
    varargout{6}=titl;
    
elseif strcmp(varargin{1},'PLOT')
    % Read the variables used to plot
    varargout{1}=xaxis;
    varargout{2}=yaxis;
    varargout{3}=titl;
    varargout{4}=bpfp;
    varargout{5}=lag_red;
    varargout{6}=datesm;
    
elseif strcmp(varargin{1},'FILTER')
    % Read the variables used to filter
    varargout{1}=xaxis;
    varargout{2}=lag_red;
    varargout{3}=fc1;
    varargout{4}=fc2;
    varargout{5}=datesm;
    varargout{6}=titl;
 
elseif strcmp(varargin{1},'DISTANCE')
    % Read the variables used to plot distance
    varargout{1}=xaxis;
    varargout{2}=vouts{3}; % pole zero file
    varargout{3}=vouts{3}; % Dateformat
    varargout{4}=lag_red;   
    varargout{5}=bpfp;  
    varargout{6}=datesm;
    
elseif strcmp(varargin{1},'CORRECT')
    % Read the variables used to CORRECT the Green's function
    varargout=vouts;
    varargout{vo+1}=filenameO;
    varargout{vo+2}=fileformatO;
    varargout{vo+3}=dateformatO;
    varargout{vo+4}=datesm;
    
end
end

