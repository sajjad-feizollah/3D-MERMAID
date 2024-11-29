% -_-_-_-_-_-_-_-_-_-_- imaging_scan_generate_Hanning -_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% generates semi-Hanning window for imaging scans.
%
% Inputs:   params:scan parameters.
% ------
% 
% Outputs:  semi-Hanning window.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function filtHann = imaging_scan_generate_Hanning(params)

HannSize=round(params.Ncol*params.HannPer)+rem(round(params.Ncol*params.HannPer),2);
filtHann=repmat(hanning(HannSize),[1,params.Npe]);
filtHann=cat(1,filtHann(1:HannSize/2,:),ones(params.Ncol-HannSize-1,params.Npe),filtHann(HannSize/2:end,:));

end