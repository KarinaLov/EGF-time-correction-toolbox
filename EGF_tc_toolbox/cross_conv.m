function [cross,varargout] = cross_conv(seq1,seq2,Fq,varargin)
  % Cross correlation function using flip and conv (convolution)
  %
  % Input: 
  %     seq1,se2 = the input segments to be cross correlated
  %     Fq = sampling frequency of the input signals
  %     varargin:
  %          wl = window length to cross correlate over in hours 
  %          swl = number of hour to stack the cross correlations over 
  %
  % Output:
  %     cross = the cross correlation of seg1 and seg 2
  %     varargout:
  %         lag = time lapse of the cross correlation
  %         faulty = vector or array giving the days were the mean are
  %                  NaN (not a number)
  %
  % Written by Karina Løviknes
  %
  
  m = length(seq1);
  n = length(seq2);
  faulty = 0;
  
  if isempty(varargin)
    % Default: the input segments are cross correlated over the entire day 
    % using flip and conv  
    s1 = flip(seq1);
    s2 = seq2;
    cross = conv(s1,s2);
    
    % Make the lag vector start from negative and have length (m+n)-1:
    k = max(m,n);
    lag = zeros(1,(m+n)-1);
    for i = 1:(m+n)-1
        lag(i) = -k+i;
    end
    
  else
    wl = varargin{1};
    swl = varargin{2};
    % The input segments are cut to shorter time windows before cross 
    % correlation to save time
    
    nwl = wl*60*60*Fq; % Window length in samples 
    
    % Preallocate for speed:
    cross_seg = zeros(swl/wl,2*nwl-1); 
    cross = zeros(24/swl,2*nwl-1);
    
    i=0;
    for l = 1:24/swl
    for j = 1:swl/wl
        i = i+1;
        s1 = seq1(1,((nwl*(i-1))+1):(nwl*i));
        s2 = seq2(1,((nwl*(i-1))+1):(nwl*i));

        sf1 = flip(s1);
        cross_seg(j,:) = conv(sf1,s2);
    end
        cross_seg_m = cross_seg(2:end-1,:);
        cross(l,:) = sum(cross_seg_m);
          % If the mean of the cross correlation is NaN the array is set to zero to avoid afffecting the final stack:
          if isnan(mean(cross(l,:))) == 1
              cross(l,:) = zeros(1,length(cross)); 
              faulty = 1;
          end
    end
    
    % Make the lag vector start from negative and have length (m+n)-1:
    mm = length(s1);
    nn = length(s2);
    k = max(mm,nn);
    lag = zeros(1,(mm+nn)-1);
    for i = 1:(mm+nn)-1
        lag(i) = -k+i;
    end
      
  end 
  
  % Outputs:
  if nargout == 2
      varargout{1} = lag;
  elseif nargout == 3
      varargout{1} = lag;
      varargout{2} = faulty;
  end
  
end