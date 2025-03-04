% -_-_-_-_-_-_-_-_-_-_-_-_-_-_- 3D_mermaid_sim -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
%
% Description: 
% -----------
% 
% Bloch simulations of 3D SE and 3D MERMAID sequences. Generates plots in
% Figures 1B, 1C, and 3.
% Each section should be independantly run.
%
% Inputs:   -
% ------
% 
% Outputs:  -
% -------
%       
% Article:  Feizollah and Tardif (2025), https://doi.org/10.1002/mrm.30436
% -------
% 
% Sajjad Feizollah, November 2024
% -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


%%  >>>>>>>>>> Bloch simulations of 3D SE and 3D MERMAID sequence <<<<<<<<<<
clear
color={
    [166,54,3]/256
    [230,85,13]/256
    [253,141,60]/256
    [0,109,44]/256
    [49,163,84]/256
    [116,196,118]/256
    [8,81,156]/256
    [49,130,189]/256
    [121, 214, 253]/256
    };
nTR=50;
TE=64;
TR=150;

% at 3T
T1=860;
T2=71;

%  >>>>>>>>>> simulate 3D SE sequence <<<<<<<<<<
label_legend{1} = "M_x";
label_legend{2} = "M_y";
label_legend{3} = "M_z";

M=bloch_simulator('se',90,180,TE,TR,nTR,T1,T2);

figure('DefaultAxesFontSize',42)
title("3D SE sequence at steady state")
hold on
plot(0:TR,M(end-TR:end,1),'LineWidth',10,'color',color{2})
plot(0:TR,M(end-TR:end,2),'LineWidth',10,'color',color{5})
plot(0:TR,M(end-TR:end,3),'LineWidth',10,'color',color{8})
legend(label_legend)
grid on
grid minor
xlabel('time [ms]','FontWeight','bold')
ylabel('magnetization','FontWeight','bold')
axis([0 TR -0.33 0.33])
axis square

%  >>>>>>>>>> simulate 3D MERMAID sequence <<<<<<<<<<
M=bloch_simulator('me',30,180,TE,TR,nTR,T1,T2,180);

figure('DefaultAxesFontSize',42)
title("3D MERMAID sequence at steady state")
hold on
plot(0:TR,M(end-TR:end,1),'LineWidth',10,'color',color{2})
plot(0:TR,M(end-TR:end,2),'LineWidth',10,'color',color{5})
plot(0:TR,M(end-TR:end,3),'LineWidth',10,'color',color{8})
legend(label_legend)
grid on
grid minor
xlabel('time [ms]','FontWeight','bold')
ylabel('magnetization','FontWeight','bold')
axis([0 TR -0.33 0.33])
axis square

%%  >>>>>>>>>> Bloch simulations for different b-values <<<<<<<<<<
clear
color={
    [166,54,3]/256
    [230,85,13]/256
    [253,141,60]/256
    [0,109,44]/256
    [49,163,84]/256
    [116,196,118]/256
    [8,81,156]/256
    [49,130,189]/256
    [121, 214, 253]/256
    };

% at 3T
T1=860;
T2=71;

gamma=267.522e6;    % gyromagnetic ratio in (rad/s/T)
G_max=73e-3; % maximum gradient amplitude for diffusion encoding (T/m)
RF_refoc_dur=5.12;  % duration of refocusing pulse
RF_exc_dur=2.56;    % duration of excitation pulse
pre_exc_dur=5.12 + 0.9 + 11 + RF_exc_dur/2;    % pre excitation duration, including inversion time, spoiler time, fat saturation, and half of excitation pulse
BW=900;
PF=6/8;
R=3;
fov=240e-3;
nTR=50;
res=(0.8:0.05:1.5)*1e-3;

for bvalue=[1000,2000,3000]*1e6
    
    for n=1:length(res)
        
        N=ceil(fov/res(n));
        [~,epi_time,dur_to_TE]=EPI_generator(fov,N,R,PF,N*BW*res(n)/res(1));
        C=[2/3*gamma^2*G_max^2,gamma^2*G_max^2*(RF_refoc_dur*1e-3+dur_to_TE*1e-3),0,-bvalue];
        root=roots(C);
        pre_readout=root(find(real(root)>0))*1000+RF_refoc_dur/2+ceil(RF_exc_dur/2);
        pre_readout_new=ceil(pre_readout)-(dur_to_TE-floor(dur_to_TE));
        TE=2*(dur_to_TE+pre_readout_new);
        TR=ceil(pre_exc_dur + TE + epi_time(end) - dur_to_TE);
        
        [~,~,Mx_se(n,bvalue/1e9)]=bloch_simulator('se',90,180,TE,TR,nTR,T1,T2);
        [~,~,Mx_pse(n,bvalue/1e9)]=bloch_simulator('me',[],180,TE,TR,nTR,T1,T2,180);
        
    end
    
end

label_legend{1} = "b = 1000 s/mm^2";
label_legend{2} = "b = 2000 s/mm^2";
label_legend{3} = "b = 3000 s/mm^2";

figure('DefaultAxesFontSize',42)
title("relative steady state signal")
hold on
plot(res*1000,Mx_pse(:,1)./Mx_se(:,1),'LineWidth',10,'color',color{2})
plot(res*1000,Mx_pse(:,2)./Mx_se(:,2),'LineWidth',10,'color',color{5})
plot(res*1000,Mx_pse(:,3)./Mx_se(:,3),'LineWidth',10,'color',color{8})
legend(label_legend)
grid on
grid minor
xlabel('resolution [mm]','FontWeight','bold')
ylabel('3D MERMAID / 3D SE','FontWeight','bold')
axis([0.7 1.6 1.5 2])
yticks(1.5:0.05:2)
xticks(0.7:0.1:1.6)
axis square


%%  >>>>>>>>>> Bloch simulations for B1+ nonuniformity <<<<<<<<<<
clear
color={
    [166,54,3]/256
    [230,85,13]/256
    [253,141,60]/256
    [0,109,44]/256
    [49,163,84]/256
    [116,196,118]/256
    [8,81,156]/256
    [49,130,189]/256
    [121, 214, 253]/256
    };

% at 3T
T1=860;
T2=71;

B1_eff=0.4:0.01:1.4;
TE=64;
TR=150;
nTR=50;

for n=1:length(B1_eff)
    [~,~,Mx_se(n)]=bloch_simulator('se',90*B1_eff(n),180*B1_eff(n),TE,TR,nTR,T1,T2);
    [~,~,Mx_pse(n)]=bloch_simulator('me',33.3*B1_eff(n),180*B1_eff(n),TE,TR,nTR,T1,T2,180*B1_eff(n));
end

label_legend{1} = "3D MERMAID";
label_legend{2} = "3D SE";

figure('DefaultAxesFontSize',42)
title("sensitivity to B1^+ nonuniformity")
hold on
plot(B1_eff,Mx_pse./max(Mx_pse),'LineWidth',10,'color',color{2})
plot(B1_eff,Mx_se./max(Mx_se),'LineWidth',10,'color',color{5})
grid on
grid minor
xlabel('relative B1^+','FontWeight','bold')
ylabel('normalized steady state signal','FontWeight','bold')
axis([0.4 1.4 0.1 1.1])
axis square
legend(label_legend)

%%  >>>>>>>>>> Bloch simulations for varying TR <<<<<<<<<<
clear
color={
    [166,54,3]/256
    [230,85,13]/256
    [253,141,60]/256
    [0,109,44]/256
    [49,163,84]/256
    [116,196,118]/256
    [8,81,156]/256
    [49,130,189]/256
    [121, 214, 253]/256
    };

% at 3T
T1wm=860;
T2wm=71;
T1gm=1300;
T2gm=72;
T1csf=4160;
T2csf=1700;

M_ss_TE = @(T1,T2,TE,TR,TI,FA)((exp(TI./T1)*(cosd(FA)-1)-exp(TR./T1)+2*exp((TE/2+TI)./T1)-cosd(FA)))...
./(-exp(TR./T1)+cosd(FA)).*exp(-TE./T2).*sind(FA); % equation (2)

TI=4;
TE=64;
TR=100:300;
FA=acosd(exp(-TR./mean([T1wm,T1gm,T1csf])));
Mx_wm=M_ss_TE(T1wm,T2wm,TE,TR,TI,FA);
Mx_gm=M_ss_TE(T1gm,T2gm,TE,TR,TI,FA);
Mx_csf=M_ss_TE(T1csf,T2csf,TE,TR,TI,FA);

label_legend{1} = "WM";
label_legend{2} = "GM";
label_legend{3} = "CSF";

figure('DefaultAxesFontSize',42)
title("3D MERMAID steady state signal vs. TR")
hold on
plot(TR,Mx_wm,'LineWidth',10,'color',color{2})
plot(TR,Mx_gm,'LineWidth',10,'color',color{5})
plot(TR,Mx_csf,'LineWidth',10,'color',color{8})
grid on
grid minor
xlabel('TR [ms]','FontWeight','bold')
ylabel('steady state signal at TE','FontWeight','bold')
axis([TR(1)-10 TR(end)+10 0 max(Mx_csf(:))+0.01])
axis square
legend(label_legend)


%%  >>>>>>>>>> Bloch simulations for varying flip angle <<<<<<<<<<
clear
color={
    [166,54,3]/256
    [230,85,13]/256
    [253,141,60]/256
    [0,109,44]/256
    [49,163,84]/256
    [116,196,118]/256
    [8,81,156]/256
    [49,130,189]/256
    [121, 214, 253]/256
    };

% at 3T
T1wm = 860;
T2wm = 71;
T1gm = 1300;
T2gm = 72;
T1csf = 4160;
T2csf =1700;

TI = 4;
TE = 64;
TR = 150;
FA = 1:0.1:90;

M_ss_TE = @(T1,T2,TE,TR,TI,FA)((exp(TI./T1)*(cosd(FA)-1)-exp(TR./T1)+2*exp((TE/2+TI)./T1)-cosd(FA)))...
./(-exp(TR./T1)+cosd(FA)).*exp(-TE./T2).*sind(FA); % equation (2)

Mx_wm=M_ss_TE(T1wm,T2wm,TE,TR,TI,FA);
Mx_gm=M_ss_TE(T1gm,T2gm,TE,TR,TI,FA);
Mx_csf=M_ss_TE(T1csf,T2csf,TE,TR,TI,FA);

label_legend{1} = "WM";
label_legend{2} = "GM";
label_legend{3} = "CSF";

figure('DefaultAxesFontSize',42)
title("steady state signal vs. flip angle")
hold on
plot(FA,Mx_wm,'LineWidth',10,'color',color{2})
plot(FA,Mx_gm,'LineWidth',10,'color',color{5})
plot(FA,Mx_csf,'LineWidth',10,'color',color{8})
grid on
grid minor
xlabel('flip angle [degree]','FontWeight','bold')
ylabel('steady state signal at TE','FontWeight','bold')
axis([0 90 0 max(Mx_csf)+0.01])
xticks(0:10:90)
axis square
legend(label_legend);
