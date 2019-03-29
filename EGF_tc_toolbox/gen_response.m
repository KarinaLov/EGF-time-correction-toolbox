function resp = gen_response(npts,sps,pz_file)
% Generates the instrument response from a pole-zero file
%
% Input:
%       npts = number of points
%       sps = samples per seconds
%       pz_file = pole zero file
%
% Output: 
%       resp = response file
%
%
% Sub-function: read_info_pz.m
%
% Written by Karina Løviknes 
% 
[pp,zz,constant] = read_info_pz(pz_file,'PZ');

% Preallocate for speed:
resp = zeros(1,npts);

for j=1:npts
    freq=j/(1/sps*npts);
    %Initialize
    freq = freq * 2 * pi;
    omega = 0 + freq * 1i;
    denom = 1+1*1i;
    num = 1+1*1i;
    
    % Zeros
    for k=1:length(zz)
        temp = omega-zz(k);
        num=num*temp;
    end
    
    % Poles
    for k=1:length(pp)
        temp = omega-pp(k);
        denom = denom * temp;
    end
    
    % Constant
    temp = real(denom) - imag(denom)*1i;
    temp = temp * num;
    mod_squared = real(denom)*real(denom)+imag(denom)*imag(denom);
    temp = real(temp) / mod_squared + (imag(temp) / mod_squared) * 1i;
    resp(j) = constant * real(temp) + constant * imag(temp) * 1i;
end
end