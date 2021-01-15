%% Calculate different BROADBAND Metrics for kronecker fixed /+ adaptive beamforming
%% Varying iterations Varying M2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;
stats_dir = 'Statistics_true_synthetic_1a' ; %%
filters_dir = 'Filters_synthetic_1a' ; %%

M = 2^6 ; 
m2_range = [1:(log2(M)-1)] ;

n_iter_range = [0, 1:2:11] ;

iSIR_dB = 0 ;

if sign( iSIR_dB ) == -1
    tmp = 'neg';
else
    tmp = '' ;
end

% Dataset
%****************************************************************************
varname = [stats_dir,'/SOI'] ;
%-----------------------------------------------------
SOI = load(varname) ; SOI = SOI.SOI ; 
num_snapshots = size(SOI,2) ; % number of speech frames
idx_snapshots = 100 ; num_snapshots ; %% snapshots = [1, 10:10:90, 100:50:450, 500:100:num_snapshots, num_snapshots] ;

num_filters = 5 ; % DS_F, MVDR_F, MVDR_K_F, DS_MVDR_F, MVDR_DS_F

intf_type = {'white', 'babble', 'hfchannel'} ;
for idx_intf_type = 1:length(intf_type)     

    Gain_dB_values = zeros( length(n_iter_range), length(m2_range), num_filters ) ;
    J_dB_values =  zeros( length(n_iter_range), length(m2_range), num_filters ) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10' ] ;
    %----------------------------------------------------- 
    % DS_F
    Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/DS_F', postname]) ; Perf = Perf.Perf ;
    Gain_dB_values(:,:,1) = Perf.Gain_dB ;
    J_dB_values(:,:,1) =  Perf.J_dB ;

    postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots) ] ;
    %----------------------------------------------------- 
    % MVDR_F
    Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_F', postname]) ; Perf = Perf.Perf ;
    Gain_dB_values(:,:,2) = Perf.Gain_dB ;
    J_dB_values(:,:,2) =  Perf.J_dB ;

    for idx_m2 = 1 : length(m2_range)

        m2 = m2_range(idx_m2) ; M2 = 2^m2 ;
        m1 = log2(M)-m2 ; M1 = 2^m1 ;

        for idx_iter = 1 : length(n_iter_range)

            n_iter = n_iter_range( idx_iter ) ;

            postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
            %-----------------------------------------------------

            % MVDR_K_F
            Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_K_F', postname]) ; Perf = Perf.Perf ;
            Gain_dB_values(idx_iter, idx_m2, 3) = Perf.Gain_dB ;
            J_dB_values(idx_iter, idx_m2, 3) =  Perf.J_dB ;

            % DS_MVDR_F
            Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/DS_MVDR_F', postname]) ; Perf = Perf.Perf ;
            Gain_dB_values(idx_iter, idx_m2, 4) = Perf.Gain_dB ;
            J_dB_values(idx_iter, idx_m2, 4) =  Perf.J_dB ;

            % MVDR_DS_F
            Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_DS_F', postname]) ; Perf = Perf.Perf ;
            Gain_dB_values(idx_iter, idx_m2, 5) = Perf.Gain_dB ;
            J_dB_values(idx_iter, idx_m2, 5) =  Perf.J_dB ;


        end %niter

    end %M2

%     % Figures
%     %-------------------------------
%     % 3 separate figures for each interference-type
%     % Each column is for a value of m2
%     figure(); 
%     subplot1( 2 , 5 , 'Min' , [0.08 0.15] , 'Max' , [0.90 1.10] , 'Gap' , [ 0.08 , 0.20] , 'XTickL' , 'All' , 'YTickL' , 'All' ) ; 
% 
%     for idx_m2 = 1 : length(m2_range)
% 
%         values = Gain_dB_values(:,idx_m2,:) ; 
%         subplot1(idx_m2) 
% 
%     %     plot(n_iter_range, values(:,:,1),'-+g'); hold on ;
%     %     plot(n_iter_range, values(:,:,2),'-xm'); hold on ;
%         plot(n_iter_range, values(:,:,3),':sr'); hold on ;
%         plot(n_iter_range, values(:,:,4),'--r'); hold on ;
%         plot(n_iter_range, values(:,:,5),'-or'); hold on ;
% 
%         ylabel('Value (dB)') ; title(['G $({\bf{h}})$ ', intf_type{idx_intf_type}]) ; axis('tight') ;
% 
%         values = J_dB_values(:,idx_m2,:) ; 
%         subplot1(idx_m2 + length(m2_range)) 
% 
%     %     plot(n_iter_range, values(:,:,1),'-+g'); hold on ;
%     %     plot(n_iter_range, values(:,:,2),'-xm'); hold on ; 
%         plot(n_iter_range, values(:,:,3),':sr'); hold on ;
%         plot(n_iter_range, values(:,:,4),'--r'); hold on ;
%         plot(n_iter_range, values(:,:,5),'-or'); hold on ;
% 
%         xlabel('$N$') ; title(['$J ({\bf{h}})$ ', intf_type{idx_intf_type} ]) ; axis('tight') ; 
% 
%     end
% 
%     hold off
% 
%     box on;
%     %legend('DS', 'MVDR', 'KP-MVDR', 'KP-DS-MVDR', 'KP-MVDR-DS') ;
%     legend('KP-MVDR', 'KP-DS-MVDR', 'KP-MVDR-DS') ;
% 
%     a=findobj(gcf); % get the handles associated with the current figure
%     allaxes=findall(a,'Type','axes');
%     alllines=findall(a,'Type','line');
%     alltext=findall(a,'Type','text');
% 
%     set(allaxes,'FontName','Times New Roman','FontWeight','Normal','LineWidth',0.5,'FontSize',16);
%     set(alllines,'Linewidth',2, 'MarkerSize', 10);
%     set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
% 
%     linkaxes( allaxes , 'xy' ) ;


    % Figures
    %-------------------------------
    % 1 figure, with each column for a type of interference
    figure(length(intf_type)+1) ;
    if idx_intf_type == 1
        subplot1( 2 , 3 , 'Min' , [0.06 0.10] , 'Max' , [1.01 1.03] , 'Gap' , [ 0.06 , 0.10] , 'XTickL' , 'All' , 'YTickL' , 'All' ) ; 
    end

    values = Gain_dB_values ; 

    subplot1(idx_intf_type) 
%     plot(m2_range, values(1,:,1),'-+g'); hold on ;
%     plot(m2_range, values(1,:,2),'-xm'); hold on ;
    plot(m2_range, values(4,:,3),':sr'); hold on ; % 4 : n_iter = 5
    plot(m2_range, values(2,:,4),'--r'); hold on ; % 2 : n_iter = 1
    plot(m2_range, values(2,:,5),'-or'); hold on ; % 2 : n_iter = 1
    title(intf_type{idx_intf_type}) ; axis('tight') ;
    ylabel('G $({\bf{h}})$ (dB)') ; 
    xticks([1:5]) ; xticklabels({'1','2','3','4','5'}) ;
%     yticks([13:18]) ; yticklabels({'13','14','15','16','17','18'}) ;

    values = J_dB_values ; 

    subplot1(idx_intf_type + length(intf_type)) 
%     plot(m2_range, values(1,:,1),'-+g'); hold on ;
%     plot(m2_range, values(1,:,2),'-xm'); hold on ; 
    plot(m2_range, values(4,:,3),':sr'); hold on ; % 4 : n_iter = 5
    plot(m2_range, values(2,:,4),'--r'); hold on ; % 2 : n_iter = 1
    plot(m2_range, values(2,:,5),'-or'); hold on ; % 2 : n_iter = 1
    ylabel('$J ({\bf{h}})$ (dB)') ; axis('tight') ;
    xlabel('$\log_2 (M_2)$') ; 
    xticks([1:5]) ; xticklabels({'1','2','3','4','5'}) ;
%     yticks([-4.5:0.5]) ; yticklabels({'-4.5','-3.5','-2.5','-1.5','-0.5','0.5'}) ;

    hold off

    if idx_intf_type == length(intf_type)
        box on;
%         legend('DS', 'MVDR', 'KP-MVDR', 'KP-DS-MVDR', 'KP-MVDR-DS') ; legend boxoff ;
        legend('KP-MVDR', 'KP-DS-MVDR', 'KP-MVDR-DS') ; legend boxoff ;

        a=findobj(gcf); % get the handles associated with the current figure
        allaxes=findall(a,'Type','axes');
        alllines=findall(a,'Type','line');
        alltext=findall(a,'Type','text');

        set(allaxes,'FontName','Times New Roman','FontWeight','Normal','LineWidth',0.5,'FontSize',16);
        set(alllines,'Linewidth',2, 'MarkerSize', 10);
        set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

        linkaxes( allaxes , 'x' ) ;
        
    end
end