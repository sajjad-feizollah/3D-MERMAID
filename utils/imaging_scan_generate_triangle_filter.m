% -_-_-_-_-_-_-_- imaging_scan_generate_triangle_filter -_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% generates triangle trajectory for imaging scan to filter the k-space for 
% subtracting phase for phase correction due to motion.
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

function filt = imaging_scan_generate_triangle_filter(params)
filt=zeros(params.Ncol,1);
filt(round(params.Ncol*abs(1-params.Rtri)/2)+1:round(params.Ncol*(1+params.Rtri)/2))=...
    [0:2/(params.Rtri*params.Ncol):1,1-2/(params.Rtri*params.Ncol):-2/(params.Rtri*params.Ncol):1/(params.Rtri*params.Ncol)];
filt=filt.*[0:2/(params.Rtri*params.Npe):1,1-2/(params.Rtri*params.Npe):-2/(params.Rtri*params.Npe):2/(params.Rtri*params.Npe)];
Ntri=size(filt,2);
filt=padarray(filt,[0,round((params.Npe-Ntri)/2)],'pre');
filt=padarray(filt,[0,floor((params.Npe-Ntri)/2)],'post');
filt(filt==0)=eps;
end