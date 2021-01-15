%% Calculate different BROADBAND Metrics for kronecker fixed /+ adaptive beamforming
%% Varying iterations Varying M2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;
filters_dir = 'Filters_synthetic_1a' ; %%

M = 2^6 ; 

intf_type = {'white', 'babble', 'hfchannel'} ;
for idx_intf_type = 1:length(intf_type) 
    
    % Dataset
    %****************************************************************************
    varname = ['Data_Snapshots/', intf_type{idx_intf_type},'/Signals_synthetic_1a'] ; %%
    %-----------------------------------------------------
    Signals = load(varname) ; 
    SOI = Signals.Signals(:,1,:) ; SOI = reshape( SOI, size(SOI,1), size(SOI,3) ) ;
    clear Signals ;
    num_snapshots = size(SOI,2) ; % number of speech frames
    snapshots = [1, 10:10:90, 100:50:450, 500:100:num_snapshots, num_snapshots] ;
    clear SOI ;

    iSIR_dB = 0 ;

    if sign( iSIR_dB ) == -1
        tmp = 'neg';
    else
        tmp = '' ;
    end

    num_filters = 5 ; % DS_F, MVDR_F, MVDR_K_F, DS_MVDR_F, MVDR_DS_F

    Gain_dB_values = zeros( length(snapshots), num_filters ) ;
    J_dB_values =  zeros( length(snapshots), num_filters ) ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10' ] ;
    %----------------------------------------------------- 
    % DS_F
    Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/DS_F', postname]) ; Perf = Perf.Perf ;
    Gain_dB_values(:,1) = Perf.Gain_dB ;
    J_dB_values(:,1) =  Perf.J_dB ;

    for idx_snaps = 1 : length(snapshots)

        idx_snapshots = snapshots(idx_snaps) ;

        postname = ['_M_', num2str(M) , '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots) ] ;
        %----------------------------------------------------- 
        % MVDR_F
        Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_F', postname]) ; Perf = Perf.Perf ;
        Gain_dB_values(idx_snaps,2) = Perf.Gain_dB ;
        J_dB_values(idx_snaps,2) =  Perf.J_dB ;

        m2 = 3 ; M2 = 2^m2 ; m1 = log2(M)-m2 ; M1 = 2^m1 ; n_iter = 5 ;
        postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
        %-----------------------------------------------------
        % MVDR_K_F
        Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_K_F', postname]) ; Perf = Perf.Perf ;
        Gain_dB_values(idx_snaps,3) = Perf.Gain_dB ;
        J_dB_values(idx_snaps,3) =  Perf.J_dB ;

        m2 = 3 ; M2 = 2^m2 ; m1 = log2(M)-m2 ; M1 = 2^m1 ; n_iter = 1 ;
        postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
        %-----------------------------------------------------
        % DS_MVDR_F
        Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/DS_MVDR_F', postname]) ; Perf = Perf.Perf ;
        Gain_dB_values(idx_snaps,4) = Perf.Gain_dB ;
        J_dB_values(idx_snaps,4) =  Perf.J_dB ;

        m2 = 3 ; M2 = 2^m2 ; m1 = log2(M)-m2 ; M1 = 2^m1 ; n_iter = 1 ;
        postname = ['_M_', num2str(M) , '_M1_' , num2str(M1) , '_M2_' , num2str(M2), '_iSIR_dB_', tmp, num2str(abs(iSIR_dB)), '_iSNR_dB_10_', 'snapshots_', num2str(idx_snapshots), '_iterations_', num2str(n_iter)] ;
        %-----------------------------------------------------
        % MVDR_DS_F
        Perf = load(['./', filters_dir,'/', intf_type{idx_intf_type}, '/MVDR_DS_F', postname]) ; Perf = Perf.Perf ;
        Gain_dB_values(idx_snaps,5) = Perf.Gain_dB ;
        J_dB_values(idx_snaps,5) =  Perf.J_dB ;

    end 


    % Figures
    %-------------------------------
    % 1 figure, with each column for a type of interference
    if idx_intf_type == 1
        figure(); %1
        subplot1( 2 , 3 , 'Min' , [0.06 0.10] , 'Max' , [1.01 1.03] , 'Gap' , [ 0.06 , 0.10] , 'XTickL' , 'All' , 'YTickL' , 'All' ) ;
    end
    
    values = Gain_dB_values ; 
    subplot1(idx_intf_type) 
    plot(snapshots, values(:,1),'-+g'); hold on ;
    plot(snapshots, values(:,2),'-xm'); 
    plot(snapshots, values(:,3),':sr');
    plot(snapshots, values(:,4),'--r');
    plot(snapshots, values(:,5),'-or');

    ylabel('G $({\bf{h}})$ (dB)') ; title(intf_type{idx_intf_type}) ; axis('tight') ;

    values = J_dB_values ; 
    subplot1(idx_intf_type + length(intf_type)) 
    plot(snapshots, values(:,1),'-+g'); hold on ;
    plot(snapshots, values(:,2),'-xm'); 
    plot(snapshots, values(:,3),':sr');
    plot(snapshots, values(:,4),'--r');
    plot(snapshots, values(:,5),'-or');

    xlabel('$r$ (snapshots)') ; ylabel('$J ({\bf{h}})$ (dB)') ; axis('tight') ;
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