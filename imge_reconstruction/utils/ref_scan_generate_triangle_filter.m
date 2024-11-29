% _-_-_-_-_-_-_-_-_- ref_scan_generate_triangle_filter _-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% generates triangle filter for phase correction.
%
% Inputs:   params: scan parameters.
% ------
% 
% Outputs:  filt: triangle filter.
% -------
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function filt = ref_scan_generate_triangle_filter(params)
filt=zeros(params.Ncol,1);

filt(round(params.Ncol*abs(1-params.Rtri)/2)+1:round(params.Ncol*(1+params.Rtri)/2))=...
    [0:2/(params.Rtri*params.Ncol):1,1-2/(params.Rtri*params.Ncol):-2/(params.Rtri*params.Ncol):1/(params.Rtri*params.Ncol)];

filt=filt.*[0:2/(params.Rtri*params.NpeRef):1,1-2/(params.Rtri*params.NpeRef):-2/(params.Rtri*params.NpeRef):2/(params.Rtri*params.NpeRef)];

filt=padarray(filt,[0,(params.NpeRef-size(filt,2))/2],'both');
end