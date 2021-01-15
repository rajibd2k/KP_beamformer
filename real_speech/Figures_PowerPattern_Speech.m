%%% Powerpatterns for kronecker fixed /+ adaptive beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;

intf_type = {'white', 'babble', 'hfchannel'} ;
idx_intf_type = 1 ;
idx_dataset = 1 ;
data_dir = ['Data/', intf_type{idx_intf_type}] ; %%
filters_dir = ['Filters_', intf_type{idx_intf_type}, '/', num2str(idx_dataset)] ; %%

iSIR_dB = 0 ; 
        
if sign( iSIR_dB ) == -1
    tmp = 'neg';
else
    tmp = '' ;
end

postname = ['_OutputVADbasedEstCurrentDisturbanceStatistics_iSNR_dB_10_iSIR_dB_', tmp, num2str(abs(iSIR_dB))] ;
%----------------------------------------------------------------------------

phi = [0:1:180]' ;
phi_rad = pi* phi / 180 ;

num_filters = 5 ;
filters = {'DS_F', 'MVDR_F', 'MVDR_K_F', 'DS_MVDR_F', 'MVDR_DS_F'} ;

num_f = 256 ;
FS = 8000 ;
f = [0 : (num_f-1)]' / num_f ; % digital frequencies 
f = f * FS ;

f_plot = [200 , 2000] ; % beampatterns will be plotted for these two frequencies - 200 Hz and 2000 Hz
[~, idx_f] = min( abs(f*ones(1,2) - ones(num_f,1)*f_plot) ) ;

PP_norm_values = zeros( length(phi), 2, num_filters ) ;

% Compute Powerpatterns of filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx_filter = 1 : num_filters
    
    filter = filters{idx_filter} ;
 
    PP = load([filters_dir,'/', filter , postname]) ; 
    
    if idx_filter == 1
        PP = PP.beta_f_1(:,idx_f) ;
    elseif idx_filter == 2
        PP = PP.beta_f_2(:,idx_f) ;
    elseif idx_filter == 3
        PP = PP.beta_f_3(:,idx_f) ;
    elseif idx_filter == 4
        PP = PP.beta_f_4(:,idx_f) ;
    elseif idx_filter == 5
        PP = PP.beta_f_5(:,idx_f) ;
    end
    
    PP = PP ./ ( ones(size(PP,1), 1) * max(PP) ) ;
    PP_norm_values(:,:,idx_filter) = PP ;
    
end

% Figures
%-------------------------------
figure(1); %1

subplot( 2,3,1 ) ; 
% DS_F
polarplot( phi_rad, PP_norm_values(:,1,1) , ':') ; grid on ; axis tight ; hold on ;
polarplot( phi_rad, PP_norm_values(:,2,1) ) ; grid on ; axis tight ; hold on ;
title('DS','FontWeight','Normal') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'','0.5',''}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
set(b,'ThetaZeroLocation' , 'Right'); 

subplot( 2,3,3 ) ; 
% MVDR_F
polarplot( phi_rad, PP_norm_values(:,1,2) , ':') ; grid on ; axis tight ; hold on ;
polarplot( phi_rad, PP_norm_values(:,2,2) ) ; grid on ; axis tight ; hold on ;
title('MVDR','FontWeight','Normal') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'','0.5',''}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
set(b,'ThetaZeroLocation' , 'Right'); 

subplot( 2,3,4 ) ; 
% MVDR_K_F
polarplot( phi_rad, PP_norm_values(:,1,3) , ':') ; grid on ; axis tight ; hold on ;
polarplot( phi_rad, PP_norm_values(:,2,3) ) ; grid on ; axis tight ; hold on ;
title('KP-MVDR','FontWeight','Normal') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'','0.5',''}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
set(b,'ThetaZeroLocation' , 'Right'); 

subplot( 2,3,5 ) ; 
% DS_MVDR_F
polarplot( phi_rad, PP_norm_values(:,1,4) , ':') ; grid on ; axis tight ; hold on ;
polarplot( phi_rad, PP_norm_values(:,2,4) ) ; grid on ; axis tight ; hold on ;
title('KP-DS-MVDR','FontWeight','Normal') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'','0.5',''}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
set(b,'ThetaZeroLocation' , 'Right'); 

subplot( 2,3,6 ) ; 
% MVDR_DS_F
polarplot( phi_rad, PP_norm_values(:,1,5) , ':') ; grid on ; axis tight ; hold on ;
polarplot( phi_rad, PP_norm_values(:,2,5) ) ; grid on ; axis tight ; hold on ;
title('KP-MVDR-DS','FontWeight','Normal') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
set(b,'RTick',[0, 0.5, 1], 'RTickLabel', {'','0.5',''}, 'ThetaLim', [0 180], 'GridAlpha', 0.5);
set(b,'ThetaTick' , [0 20 30 70 130 160 180],'ThetaTickLabel' , {'0^o', '','','70^o', '','','180^o'}); 
set(b,'ThetaZeroLocation' , 'Right'); 

box on ;
legend('200 Hz', '2000 Hz') ; legend boxoff ;

a=findobj(gcf); % get the handles associated with the current figure

alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

