% -_-_-_-_-_-_-_-_-_-_-_-_-_-_- bloch_simulator -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% Simulates evolution of magnetizations of 3D SE and 3D MERMAID sequences
% using Bloch equations
%
% Inputs:
% ------
% 
%     type: type of sequence:
%                   - 'se' for convetional 3D spin-echo sequence
%                   - 'me' for proposed sequence
% 
%     exc: excitation flip angle
% 
%     refocusing: refocusing pulse flip angle
% 
%     TE: echo time of the sequence (ms)
% 
%     TR: repetition time (ms)
%       
%     nTR: number of TRs to simulate
%
%     inversion:  inversion flip angle for 'pse' sequence
% 
% Outputs:
% -------
% 
%    M: simulated magnetization over multiple TRs
% 
%    Mz_ss: longitudinal magnetization at steady state
% 
%    Mx_ss: transverse magnetization at steady state
%       
% Article: 
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

function [M,Mz_ss,Mx_ss]=bloch_simulator(type,exc,refoc,TE,TR,nTR,T1,T2,inversion)

% >>>>>>>>>> sequence parameters <<<<<<<<<<

dur_exc=2;         % duration of excitation pulse
dur_inv=6;         % duration of the inversion pulse
dur_spoiler=1;     % duration of spoiler before excitation
M=[0,0,1]';        % initial magnetization
if(isempty(exc))
    exc=acosd(exp(-TR/T1));
end
TE=2*floor(TE/2);  % only even TEs are used

% >>>>>>>>>> simulation parameters <<<<<<<<<<

[a,b]=freeprecess(1,T1,T2,0);  % precession function

% >>>>>>>>>> Bloch simulations <<<<<<<<<<

if(type=="me")     % proposed sequence
    
    for n=1:nTR
        
        % start of sequence to inversion
        for k=(n-1)*TR:(n-1)*TR+dur_inv/2
            M(:,k+2)=a*M(:,end)+b;
        end
        
        % inversion to excitation
        M=[M,xrot(inversion)*M(:,end)];
        for k=(n-1)*TR+dur_inv/2+1:(n-1)*TR+dur_inv/2+dur_spoiler+dur_exc/2+1
            M(:,k+2)=a*M(:,end)+b;
        end
        M(2,:)=0;
        
        % excitation to refocusing
        M=[M,yrot(exc)*M(:,end)];
        for k=(n-1)*TR+dur_inv/2+dur_spoiler+dur_exc/2+2:...
                (n-1)*TR+dur_inv/2+dur_spoiler+dur_exc/2+TE/2+2
            M(:,k+2)=a*M(:,end)+b;
        end
        
        % refocusing to end of TR
        M=[M,xrot(-refoc)*M(:,end)];
        for k=(n-1)*TR+dur_inv/2+dur_spoiler+dur_exc/2+TE/2+3:...
                n*TR-1
            M(:,k+2)=a*M(:,end)+b;
        end
        
        % spoiling of transverse magnetization
        M(2,:)=0;
        M(1:2,end)=0;
    end
    
elseif (type=="se") % conventional 3D SE sequence
    
    for n=1:nTR
        
        % start of sequence to excitation
        for k=(n-1)*TR:(n-1)*TR+dur_exc/2
            M(:,k+2)=a*M(:,end)+b;
        end
        
        % excitation to refocusing
        M=[M,yrot(exc)*M(:,end)];
        for k=(n-1)*TR+dur_exc/2+1:(n-1)*TR+dur_exc/2+TE/2+1
            M(:,k+2)=a*M(:,end)+b;
        end
        
        % refocusing to end of TR
        M=[M,xrot(refoc)*M(:,end)];
        for k=(n-1)*TR+dur_exc/2+TE/2+2:n*TR-1
            M(:,k+2)=a*M(:,end)+b;
        end
        
        % spoiling of transverse magnetization
        M(2,:)=0;
        M(1:2,end)=0;
    end
else
    error('Sequence not defined!')
end

% calculate steady-state magnetization
Mz_ss=max(M(3,(nTR-1)*TR-10:end));
Mx_ss=sind(exc)*Mz_ss*exp(-TE/T2);
M=M';
end

% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%  Functions required to calculate rotation matrices.
%  Adapted from Brian Hargreaves' online tutorial:
%  mrsrl.stanford.edu/~brian/bloch/
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

% ----------------------------------------------------------------
% calculates A and B for magnetization in M=A*M0+B for the time T
% ----------------------------------------------------------------
function [Afp,Bfp]=freeprecess(T,T1,T2,df)
M0=1;
E1=exp(-T/T1);
E2=exp(-T/T2);
Afp=[E2 0 0;0 E2 0;0 0 E1]*zrot((360*df)*(T/1000)); 
Bfp=[0 0 M0-M0*E1]';
end

% ----------------------------------------------------------------
% rotation matrix around X axis
% ----------------------------------------------------------------
function Rx=xrot(phi)
Rx = [1 0 0; 0 cosd(phi) -sind(phi);0 sind(phi) cosd(phi)];
end

% ----------------------------------------------------------------
% rotation matrix around Y axis
% ----------------------------------------------------------------
function Ry=yrot(phi)
Ry = [cosd(phi) 0 sind(phi);0 1 0;-sind(phi) 0 cosd(phi)];
end

% ----------------------------------------------------------------
% rotation matrix around Z axis
% ----------------------------------------------------------------
function Rz=zrot(phi)
Rz=[cosd(phi) -sind(phi) 0;sind(phi) cosd(phi) 0;0, 0, 1];
end
