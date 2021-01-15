%%% Calculate Powerpatterns for kronecker fixed /+ adaptive beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;
stats_dir = 'Statistics_true_synthetic_1a' ; %%
filters_dir = 'Filters_synthetic_1a' ; %%

theta_d = 70 ;
theta_u = [20, 30, 130, 160] ; 

delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
c = 340 ; % m/s
Ts = 1/8000 ; % s

M = 2^6 ;

iSIR_dB = 0 ;

if sign( iSIR_dB ) == -1
    tmp = 'neg';
else
    tmp = '' ;
end

phi=0:1:180 ;
phi_rad=phi/180*pi;

nfft = 256 ; % 256 pt fft
f = [0 : (nfft-1)]' / nfft ;

num_filters = 5 ;

PP_norm_values = zeros( length(phi_rad), num_filters ) ;

% Dataset
%****************************************************************************
varname = [stats_dir,'/SOI'] ;
%-----------------------------------------------------
SOI = load(varname) ; SOI = SOI.SOI ; 
num_snapshots = size(SOI,2) ; % number of speech frames
idx_snapshots = num_snapshots ;

% Compute Powerpatterns of filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10' ] ;
%----------------------------------------------------- 
h = load(['./', filters_dir,'/DS_F', postname]) ; h = h.h ;
MeanPPf = zeros(length(phi),length(f)) ; 
for idx_f = 1:length(f)
    [~, ~, PP_f_norm] = PP( f(idx_f), h(:,idx_f), delta, c, Ts) ;
    MeanPPf(:,idx_f) = PP_f_norm ;
end
MeanPPf = nanmean( MeanPPf, 2 ) ; MeanPPf = MeanPPf / max(MeanPPf) ; 
PP_norm_values(:,1) = MeanPPf ;



postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots) ] ;
%-----------------------------------------------------
h = load(['./', filters_dir, '/MVDR_F', postname]) ; h = h.h ;
MeanPPf = zeros(length(phi),length(f)) ; 
for idx_f = 1:length(f)
    [~, ~, PP_f_norm] = PP( f(idx_f), h(:,idx_f), delta, c, Ts) ;
    MeanPPf(:,idx_f) = PP_f_norm ;
end
MeanPPf = nanmean( MeanPPf, 2 ) ; MeanPPf = MeanPPf / max(MeanPPf) ; 
PP_norm_values(:,2) = MeanPPf ;



m2 = 3 ; M2 = 2^m2 ; M1 = M / M2 ; n_iter = 5 ;
varname = [stats_dir,'/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ;
postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
%----------------------------------------------------- 
h = load(['./', filters_dir, '/MVDR_K_F', postname]) ; h = h.h ;
MeanPPf = zeros(length(phi),length(f)) ; 
for idx_f = 1:length(f)
    [~, ~, PP_f_norm] = PP( f(idx_f), h(:,idx_f), delta, c, Ts) ;
    MeanPPf(:,idx_f) = PP_f_norm ;
end
MeanPPf = nanmean( MeanPPf, 2 ) ; MeanPPf = MeanPPf / max(MeanPPf) ; 
PP_norm_values(:,3) = MeanPPf ;



m2 = 5 ; M2 = 2^m2 ; M1 = M / M2 ; n_iter = 1 ;
varname = [stats_dir,'/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ;
postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
%----------------------------------------------------- 
h = load(['./', filters_dir, '/DS_MVDR_F', postname]) ; h = h.h ;
MeanPPf = zeros(length(phi),length(f)) ; 
for idx_f = 1:length(f)
    [~, ~, PP_f_norm] = PP( f(idx_f), h(:,idx_f), delta, c, Ts) ;
    MeanPPf(:,idx_f) = PP_f_norm ;
end
MeanPPf = nanmean( MeanPPf, 2 ) ; MeanPPf = MeanPPf / max(MeanPPf) ; 
PP_norm_values(:,4) = MeanPPf ;



m2 = 3 ; M2 = 2^m2 ; M1 = M / M2 ; n_iter = 1 ;
varname = [stats_dir,'/d_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2) ] ;
Steering_vectors = load(varname) ; d_1 = Steering_vectors.d_1 ; d_2 = Steering_vectors.d_2 ;
postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
%----------------------------------------------------- 
h = load(['./', filters_dir, '/MVDR_DS_F', postname]) ; h = h.h ;
MeanPPf = zeros(length(phi),length(f)) ; 
for idx_f = 1:length(f)
    [~, ~, PP_f_norm] = PP( f(idx_f), h(:,idx_f), delta, c, Ts) ;
    MeanPPf(:,idx_f) = PP_f_norm ;
end
MeanPPf = nanmean( MeanPPf, 2 ) ; MeanPPf = MeanPPf / max(MeanPPf) ; 
PP_norm_values(:,5) = MeanPPf ;


% Figures
%-------------------------------
figure(); %1

% subplot( 2,3,1 ) ; 
% % DS_F
% polarplot( phi_rad, PP_norm_values(:,1) ) ; grid on ; axis tight ; hold on ;
% title('DS','FontWeight','Normal') ;
% b=gca;
% set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
% set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'0','0.5','1'}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
% set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
% set(b,'ThetaZeroLocation' , 'Right'); 

subplot( 2,2,1 ) ; 
% MVDR_F
polarplot( phi_rad, PP_norm_values(:,2) ) ; grid on ; axis tight ; hold on ;
title('MVDR','FontWeight','Normal') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'0','0.5','1'}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
set(b,'ThetaZeroLocation' , 'Right'); 

subplot( 2,2,2 ) ; 
% MVDR_K_F
polarplot( phi_rad, PP_norm_values(:,3) ) ; grid on ; axis tight ; hold on ;
title('KP-MVDR','FontWeight','Normal') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'0','0.5','1'}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
set(b,'ThetaZeroLocation' , 'Right'); 

subplot( 2,2,3 ) ; 
% DS_MVDR_F
polarplot( phi_rad, PP_norm_values(:,4) ) ; grid on ; axis tight ; hold on ;
title('KP-DS-MVDR','FontWeight','Normal') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'0','0.5','1'}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
set(b,'ThetaZeroLocation' , 'Right'); 

subplot( 2,2,4 ) ; 
% MVDR_DS_F
polarplot( phi_rad, PP_norm_values(:,5) ) ; grid on ; axis tight ; hold on ;
title('KP-MVDR-DS','FontWeight','Normal') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'0','0.5','1'}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
set(b,'ThetaZeroLocation' , 'Right'); 

box on ;

a=findobj(gcf); % get the handles associated with the current figure

alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

