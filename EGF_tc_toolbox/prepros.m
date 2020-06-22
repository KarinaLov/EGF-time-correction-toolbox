function prosd = prepros(data,Fq,df,resp,channel,varargin)
% Preprocesses the data before cross correlations
%
% Input: 
%       data = input data to be processed
%       Fq = sampling frequency
%       df = designfilter
%       pz_pile = pole-zero file
%       varargin: inputs about what normalization processes to apply 
%                 default: onebit normalization before spectral whitening
%                 'onebit', 'sw': onebit normalization before spectral whitening
%                 'onebit': only onebit normalization
%                 'sw': spectral whitening before onebit normalization
%
% Output: 
%       prosd = the processed data
%
%
% Sub-function: costap_filter.m, rm_resp.m, spectral_whitening.m
%
% Written by Karina Løviknes
%

L = length(data);

% Remove mean and trend
dtrnd = detrend(data);

% Taper
tap = costap_filter(dtrnd,0.05);

% Bandpass
filt1 = filtfilt(df,tap);
% Remove instrument response 
trans = rm_iresp(filt1,L,resp);
% Bandpass again to avoid low frequency artifacts caused by 
filt2 = filtfilt(df,trans);

if isempty(varargin)
    % Default: Onebit normalization before spectral whitening
    onebit = filt2./abs(filt2);        
    white = spectral_whitening(onebit);    
    prosd = white;
else 
    norm = varargin{1};
    
    if length(norm) == 1 && strcmp(norm{1},'onebit')
        % Only onebit normalization
        onebit = filt2./abs(filt2);        
        prosd = onebit;

    elseif length(norm) == 1 && strcmp(norm{1},'sw')
        % Spectral whitening before onebit normalize     
        white = spectral_whitening(filt2);
        onebit = white./abs(white);        
        prosd = onebit;

    elseif strcmp(norm{1},'onebit') && strcmp(norm{2},'sw')
        % Default: Onebit normalization before spectral whitening
        onebit=filt2./abs(filt2);        
        white=spectral_whitening(onebit);    
        prosd=white;

    % elseif strcmp(varargin{1},'rms')
    %        rms normalization (not yet an option) 
    end
end                    
end
