function varargout = read_info_pz(pole_zero_filename,varargin)
% Reads the commented information in the pole zero files
%
% Input:
%       pole_zero_file_name = the name of the settingsfile
%       varargin: spesfification about the type of values to read
%               'PZ': poles, zeros and constant
%               'Coordinates': latitude and longitude
%               'Station': network, stationname, location and channels
%
% Output:
%       varargout
%
% 

% Open the file
pz_fid = fopen(pole_zero_filename,'r');
%txt = readtable(pole_zero_file_name,'FileType','text','Delimiter',':','Delimiter',' ');

% Loop over the entire file:
line = 0;
while (line ~= -1)
    % Read the next line in the file:
    line = fgets(pz_fid);
    % Check to make sure it is not the end of the file:
    if (line ~= -1)

        txtcell = textscan(line,'%s','Delimiter',':');
        txt=txtcell{1};

        if length(txt) > 1
            var_name=txt{1};
            var=txt{2};
        
            if strcmp(var_name,'* NETWORK   (KNETWK)')
                network = var;
                    
            elseif strcmp(var_name,'* STATION    (KSTNM)')
                station = var;
                    
            elseif strcmp(var_name,'* LOCATION   (KHOLE)')
                location = var;
                    
            elseif strcmp(var_name,'* CHANNEL   (KCMPNM)')
                channel = var;
                    
            elseif strcmp(var_name,'* CREATED           ')
                var = [txt{2:end}];
                created = var;
                    
            elseif strcmp(var_name,'* START             ')
                var = [txt{2:end}];
                starttime = var;
                    
            elseif strcmp(var_name,'* END               ')
                var = [txt{2:end}];
                endtime = var;
                    
            elseif strcmp(var_name,'* DESCRIPTION       ')
                created = var;
                            
            elseif strcmp(var_name, '* LATITUDE    (deg) ')
                latitude = str2num(var);
                    
            elseif strcmp(var_name,'* LONGITUDE   (deg) ')
                longitude = str2num(var);
                    
            elseif strcmp(var_name,'* ELEVATION     (m) ')
                elevation = str2num(var);
                    
            elseif strcmp(var_name,'* DEPTH         (m) ')
                depth = str2num(var);
                    
            elseif strcmp(var_name,'* DIP         (deg) ')
                dip = str2num(var);
                    
            elseif strcmp(var_name,'* AZIMUTH     (deg) ')
                azimuth = str2num(var);
                    
            elseif strcmp(var_name,'* SAMPLE RATE  (Hz) ')
                samplerate = str2num(var);
                    
            elseif strcmp(var_name,'* INPUT UNIT        ')
                input_units = var;
                    
            elseif strcmp(var_name,'* OUTPUT UNIT       ')
                output_units = var;
                    
            elseif strcmp(var_name,'* INSTTYPE          ')
                insttype = var;
                    
            elseif strcmp(var_name,'* INSTGAIN          ')
                instgain = var;
                    
            elseif strcmp(var_name,'* COMMENT           ')
                comment = var;
                    
            elseif strcmp(var_name,'* SENSITIVITY       ')
                sensitivity = str2num(var);
                    
            elseif strcmp(var_name,'* A0                ')
                AO = str2num(var);
                    
            elseif strcmp(var_name,'* Site Name         ')
                site_name = var;
                    
            elseif strcmp(var_name,'* Owner             ')
                owner = var;
            end
        elseif length(txt)==1
            txt2 = split(txt);
            var_name = txt2{1};
            var = txt2{2};
            
            if strcmp(var_name,'ZEROS')
                num_zeros = str2double(var);
                zzs = zeros(1,num_zeros);
                for k = 1:num_zeros                   
                    line0 = fgets(pz_fid);

                    txtcell0 = textscan(line0,'%s');
                    txt0 = split(txtcell0{1});

                    zzs(k) = str2double(txt0{1}) + str2double(txt0{2})*1i;
                end              
            elseif strcmp(var_name,'POLES')
                num_poles = str2double(var);
                pps = zeros(1,num_poles);
                for k = 1:num_poles                   
                    linep = fgets(pz_fid);

                    txtcellp = textscan(linep,'%s');
                    txtp = split(txtcellp{1});

                    pps(k) = str2double(txtp{1}) + str2double(txtp{2})*1i;
                end          
            elseif strcmp(var_name,'CONSTANT')
                constant = str2double(var);
            end                
            end
        end
end

if strcmp(varargin{1},'PZ')
    varargout{1} = pps;
    varargout{2} = zzs;
    varargout{3} = constant;
    
elseif strcmp(varargin{1},'Coordinates')
    varargout{1} = latitude;
    varargout{2} = longitude;
        
elseif strcmp(varargin{1},'Station')
    varargout{1} = network;
    varargout{2} = station;
    varargout{3} = location;
    varargout{4} = channel;
end
fclose(pz_fid);
end
    
    