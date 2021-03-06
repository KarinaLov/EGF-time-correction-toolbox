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
            var_name = txt{1};
            var = txt{2};
        
            if strcmpi(var_name,'* NETWORK   (KNETWK)')
                network = var;
                    
            elseif strcmpi(var_name,'* STATION    (KSTNM)')
                station = var;
                    
            elseif strcmpi(var_name,'* LOCATION   (KHOLE)')
                location = var;
                    
            elseif strcmpi(var_name,'* CHANNEL   (KCMPNM)')
                channel = var;
                    
            elseif strcmpi(var_name,'* CREATED           ')
                var = [txt{2:end}];
                created = var;
                    
            elseif strcmpi(var_name,'* START             ')
                var = [txt{2:end}];
                starttime = var;
                    
            elseif strcmpi(var_name,'* END               ')
                var = [txt{2:end}];
                endtime = var;
                    
            elseif strcmpi(var_name,'* DESCRIPTION       ')
                created = var;
                            
            elseif strcmpi(var_name, '* LATITUDE    (deg) ') ||...
                    strcmp(var_name, '* Latitude          ')
                latitude = str2num(var);
                    
            elseif strcmpi(var_name, '* LONGITUDE   (deg) ') ||...
                    strcmp(var_name, '* Longitude         ')
                longitude = str2num(var);
                    
            elseif strcmpi(var_name,'* ELEVATION     (m) ') ||...
                    strcmp(var_name,'* Elevation         ') 
                elevation = str2num(var);
                    
            elseif strcmpi(var_name,'* DEPTH         (m) ') ||...
                    strcmp(var_name,'* Depth             ') 
                depth = str2num(var);
                    
            elseif strcmpi(var_name,'* DIP         (deg) ') ||...
                    strcmp(var_name,'* Dip               ')
                dip = str2num(var);
                    
            elseif strcmpi(var_name,'* AZIMUTH     (deg) ') ||...
                    strcmp(var_name,'* Azimuth           ')
                azimuth = str2num(var);
                    
            elseif strcmpi(var_name,'* SAMPLE RATE  (Hz) ') ||...
                    strcmp(var_name,'* SAMPLE RATE       ')
                samplerate = str2num(var);
                    
            elseif strcmpi(var_name,'* INPUT UNIT        ')
                input_units = var;
                    
            elseif strcmpi(var_name,'* OUTPUT UNIT       ')
                output_units = var;
                    
            elseif strcmpi(var_name,'* INSTTYPE          ')
                insttype = var;
                    
            elseif strcmpi(var_name,'* INSTGAIN          ')
                instgain = var;
                    
            elseif strcmpi(var_name,'* COMMENT           ')
                comment = var;
                    
            elseif strcmpi(var_name,'* SENSITIVITY       ')
                sensitivity = str2num(var);
                    
            elseif strcmpi(var_name,'* A0                ')
                AO = str2num(var);
                    
            elseif strcmpi(var_name,'* Site Name         ')
                site_name = var;
                    
            elseif strcmpi(var_name,'* Owner             ')
                owner = var;
            end
        elseif length(txt)==1
            txt2 = split(txt);
            var_name = txt2{1};
            var = txt2{2};
            
        if strcmpi(var_name,'ZEROS')
                num_zeros = str2double(var);
                zzs = zeros(1,num_zeros);
                pos01 = ftell(pz_fid);
                k = 1;
                while k < num_zeros              
                    line0 = fgets(pz_fid);
                    pos0(k) = ftell(pz_fid);
                    if ~strncmp(line0,'P',1) 
                        txtcell0 = textscan(line0,'%s');
                        txt0 = split(txtcell0{1});

                        if ~isnan(str2double(txt0{1}))
                            zzs(k) = str2double(txt0{1}) + str2double(txt0{2})*1i;
                            k = k+1;            
                        else                       
                            fseek(pz_fid,pos0(k-1),'bof');
                            k = num_zeros+1;
                        end
                    else
                        fseek(pz_fid,pos01,'bof');
                         k = num_zeros+1;                        
                    end
                end
            elseif strcmpi(var_name,'POLES')                                
                num_poles = str2double(var);
                pps = zeros(1,num_zeros);
                k = 1;
                while k < num_poles                     
                    linep = fgets(pz_fid);
                    posp(k) = ftell(pz_fid);

                    txtcellp = textscan(linep,'%s');
                    txtp = split(txtcellp{1});
                    
                    if ~isnan(str2double(txtp{1}))
                        pps(k) = str2double(txtp{1}) + str2double(txtp{2})*1i;
                        k = k+1;              
                    else
                        fseek(pz_fid,posp(k-1),'bof');
                        k = num_poles+1;
                    end
                end
            elseif strcmp(var_name,'CONSTANT')                                
                constant = str2double(var);
            end
        end
    end
end

if strcmpi(varargin{1},'PZ')
    varargout{1} = pps;
    varargout{2} = zzs;
    varargout{3} = constant;
    
elseif strcmpi(varargin{1},'Coordinates')
    varargout{1} = latitude;
    varargout{2} = longitude;
        
elseif strcmpi(varargin{1},'Station')
    varargout{1} = network;
    varargout{2} = station;
    varargout{3} = location;
    varargout{4} = channel;
end
fclose(pz_fid);
end
    
    