%% Misalignment Errors for kronecker fixed /+ adaptive beamforming
%% Varying iSIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;

M = 2^6 ; 

iSIR_dB_range = [-10:5:20] ;
num_filters = 5 ;
intf_type = {'white', 'babble', 'hfchannel'} ;

Perf_Metrics = cell(1, length(intf_type)) ;
for idx_intf_type = 1:length(intf_type) 
    
    % Compute MisAlignment Errors of filters, PROCESSED and UNPROCESSED output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Perf = zeros( length(iSIR_dB_range), num_filters , 2 ) ;

    for idx_iSIR_dB = 1 : length(iSIR_dB_range)

        iSIR_dB = iSIR_dB_range(idx_iSIR_dB) ;

        if sign( iSIR_dB ) == -1
            tmp = 'neg';
        else
            tmp = '' ;
        end
               
        for idx_dataset = [1:8,11:18] %% 1:100
            
            filters_dir = ['Filters_', intf_type{idx_intf_type}, '/', num2str(idx_dataset)] ; %%

            postname = ['_OutputVADbasedEstCurrentDisturbanceStatistics_iSNR_dB_10_iSIR_dB_', tmp, num2str(abs(iSIR_dB))] ;
            %----------------------------------------------------------------------------

            filter = 'DS_F' ; idx_filter = 1 ;
            Metric = load([filters_dir,'/', filter , postname]) ; 
            Metric = Metric.MLA; 
            Metric = reshape( Perf(idx_iSIR_dB,idx_filter,:), 1, 2 ) + Metric ;
            Perf(idx_iSIR_dB,idx_filter,:) = Metric ;

            filter = 'MVDR_F' ; idx_filter = 2 ;
            Metric = load([filters_dir,'/', filter , postname]) ; 
            Metric = Metric.MLA; 
            Metric = reshape( Perf(idx_iSIR_dB,idx_filter,:), 1, 2 ) + Metric ;
            Perf(idx_iSIR_dB,idx_filter,:) = Metric ;

            filter = 'MVDR_K_F' ; idx_filter = 3 ;
            Metric = load([filters_dir,'/', filter , postname]) ; 
            Metric = Metric.MLA; 
            Metric = reshape( Perf(idx_iSIR_dB,idx_filter,:), 1, 2 ) + Metric ;
            Perf(idx_iSIR_dB,idx_filter,:) = Metric ;

            filter = 'DS_MVDR_F' ; idx_filter = 4 ;
            Metric = load([filters_dir,'/', filter , postname]) ; 
            Metric = Metric.MLA; 
            Metric = reshape( Perf(idx_iSIR_dB,idx_filter,:), 1, 2 ) + Metric ;
            Perf(idx_iSIR_dB,idx_filter,:) = Metric ;

            filter = 'MVDR_DS_F' ; idx_filter = 5 ;
            Metric = load([filters_dir,'/', filter , postname]) ; 
            Metric = Metric.MLA; 
            Metric = reshape( Perf(idx_iSIR_dB,idx_filter,:), 1, 2 ) + Metric ;
            Perf(idx_iSIR_dB,idx_filter,:) = Metric ;
        
        end

    end 
    
    Perf = Perf / 16 ; % 16 datasets
    Perf_Metrics{idx_intf_type} = Perf ;

end

save('MLA', 'Perf_Metrics') ;

clear all ; clc ; close all ;
iSIR_dB_range = [-10:5:20] ;
num_filters = 5 ;
intf_type = {'white', 'babble', 'hfchannel'} ;
Perf_Metrics = load('MLA') ;
Perf_Metrics = Perf_Metrics.Perf_Metrics ;

% Figures
%-------------------------------
% 1 figure, with each column for a type of interference
figure(); %1
subplot1( 1 , 3 , 'Min' , [0.06 0.10] , 'Max' , [1.01 1.03] , 'Gap' , [ 0.06 , 0.10] , 'XTickL' , 'All' , 'YTickL' , 'All' ) ;

for idx_intf_type = 1:length(intf_type)  
    
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

    xlabel('iSIR (dB)') ; ylabel('MLA (dB)') ; title(intf_type{idx_intf_type}) ; axis('tight') ;
    xticks([-10:10:20]) ; xticklabels({'-10','0','10','20'}) ;
    %yticks([13:20]) ; yticklabels({'13','14','15','16','17','18','19','20'}) ;

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
