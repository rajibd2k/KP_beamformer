%%% Metrics (Misalignment, PESQ-BB, STOI) for kronecker fixed /+ adaptive beamforming
%%% Varying iSIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; clc ; close all ;
iSIR_dB_range = [-10:5:20] ;
num_filters = 5 ;
intf_type = {'white', 'babble', 'hfchannel'} ;

% Figures
%-------------------------------
% 1 figure, with each column for a type of interference
figure(); %1
subplot1( 3 , 3 , 'Min' , [0.05 0.09] , 'Max' , [1.03 0.99] , 'Gap' , [ 0.08 , 0.04] , 'XTickL' , 'Margin' , 'YTickL' , 'All' ) ;

for idx_intf_type = 1:length(intf_type)  
    
    % MLA
    Perf_Metrics = load('MLA') ;
    Perf_Metrics = Perf_Metrics.Perf_Metrics ;
    
    values = Perf_Metrics{idx_intf_type} ;

    % unprocessed or processed output
    values = values(:,:,1 ) ; %% 1 - Output or 2 - Output_postprocess
    values = reshape( values, length(iSIR_dB_range), num_filters ) ; 
    values = 10*log10( values/100 ) ; % MLA was calculated as percentage

    subplot1(idx_intf_type) 
    plot(iSIR_dB_range, values(:,1),'-+g'); hold on ;
    plot(iSIR_dB_range, values(:,2),'-xm'); 
    plot(iSIR_dB_range, values(:,3),':sr');
    plot(iSIR_dB_range, values(:,4),'--r');
    plot(iSIR_dB_range, values(:,5),'-or');

    ylabel('Misalignment (dB)') ; title(intf_type{idx_intf_type}) ; axis('tight') ;
    xlim([-10,10]) ; %xticks([-10:5:10]) ; xticklabels({'-10','-5', '0', '5', '10'}) ;
    hold off
    
    % PESQ-BB
    Perf_Metrics = load('PESQ_BB') ;
    Perf_Metrics = Perf_Metrics.Perf_Metrics ;
    
    values = Perf_Metrics{idx_intf_type} ;

    % unprocessed or processed output
    values = values(:,:,1 ) ; %% 1 - Output or 2 - Output_postprocess
    values = reshape( values, length(iSIR_dB_range), num_filters ) ; 
    values = 10*log10( values ) ; 

    subplot1(idx_intf_type + length(intf_type)) 
    plot(iSIR_dB_range, values(:,1),'-+g'); hold on ;
    plot(iSIR_dB_range, values(:,2),'-xm'); 
    plot(iSIR_dB_range, values(:,3),':sr');
    plot(iSIR_dB_range, values(:,4),'--r');
    plot(iSIR_dB_range, values(:,5),'-or');

    ylabel('PESQ (dB)') ; axis('tight') ;
    xlim([-10,10]) ; %xticks([-10:5:10]) ; xticklabels({'-10','-5', '0', '5', '10'}) ;
    hold off
    
    
    % STOI
    Perf_Metrics = load('STOI') ;
    Perf_Metrics = Perf_Metrics.Perf_Metrics ;
    
    values = Perf_Metrics{idx_intf_type} ;

    % unprocessed or processed output
    values = values(:,:,1 ) ; %% 1 - Output or 2 - Output_postprocess
    values = reshape( values, length(iSIR_dB_range), num_filters ) ; 
    values = 10*log10( values ) ; 

    subplot1(idx_intf_type + 2*length(intf_type)) 
    plot(iSIR_dB_range, values(:,1),'-+g'); hold on ;
    plot(iSIR_dB_range, values(:,2),'-xm'); 
    plot(iSIR_dB_range, values(:,3),':sr');
    plot(iSIR_dB_range, values(:,4),'--r');
    plot(iSIR_dB_range, values(:,5),'-or');

    xlabel('iSIR (dB)') ; ylabel('STOI (dB)') ; axis('tight') ;
    xlim([-10,10]) ; xticks([-10:5:10]) ; xticklabels({'-10','-5', '0', '5', '10'}) ;
    hold off

end

box on;
legend('DS', 'MVDR', 'KP-MVDR', 'KP-DS-MVDR', 'KP-MVDR-DS') ; legend boxoff ;

a=findobj(gcf); % get the handles associated with the current figure
allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Times New Roman','FontWeight','Normal','LineWidth',0.5,'FontSize',16);
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

linkaxes( allaxes , 'x' ) ;
