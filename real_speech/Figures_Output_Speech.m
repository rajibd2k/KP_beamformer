%% Outputs of kronecker fixed /+ adaptive beamforming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; close all ; clc ;

intf_type = {'white', 'babble', 'hfchannel'} ;
idx_intf_type = 1 ;
idx_dataset = 1 ;
data_dir = ['Data/', intf_type{idx_intf_type}] ; %%
filters_dir = ['Filters_', intf_type{idx_intf_type}, '/', num2str(idx_dataset)] ; %%

SOI = load([data_dir,'/SOI_' , num2str(idx_dataset)]) ; SOI = SOI.SOI(:,:,1) ; SOI = reshape( SOI, size(SOI,1) , 1 ) ;
FS = 8000 ; 
tms = [0 : (length(SOI) - 1)]' / FS ; % s

iSIR_dB = 0 ; 
        
if sign( iSIR_dB ) == -1
    tmp = 'neg';
else
    tmp = '' ;
end

postname = ['_', num2str(idx_dataset), '_iSNR_dB_10_iSIR_dB_', tmp, num2str(abs(iSIR_dB))] ;
%----------------------------------------------------------------------------
Input = load( [data_dir,'/Input', postname] ) ; 
Input = Input.Input(:,1) ;

postname = ['_OutputVADbasedEstCurrentDisturbanceStatistics_iSNR_dB_10_iSIR_dB_', tmp, num2str(abs(iSIR_dB))] ;
%----------------------------------------------------------------------------

figure(1); %1
subplot1( 7 , 1 , 'Min' , [0.05 0.07] , 'Max' , [0.8 1.01] , 'Gap' , [ 0.00 , 0.04] , 'XTickL' , 'Margin' , 'YTickL' , 'All' ) ; 

subplot1(1) ; plot(tms, SOI) ; axis('tight') ; title('$x(t)$','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;

subplot1(2) ; plot(tms, Input) ; axis('tight') ; title('$y_1(t)$','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;
   
filter = 'DS_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output ; % signal = signal.Output_postprocess ;
subplot1(3) ; plot(tms, signal) ; axis('tight') ; title('$z(t)$ of DS','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;

filter = 'MVDR_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output ; % signal = signal.Output_postprocess ;
subplot1(4) ; plot(tms, signal) ; axis('tight') ; title('$z(t)$ of MVDR','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;
ylabel('Amplitude') ;

filter = 'MVDR_K_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output ; % signal = signal.Output_postprocess ;
subplot1(5) ; plot(tms, signal) ; axis('tight') ; title('$z(t)$ of KP-MVDR','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;

filter = 'DS_MVDR_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output ; % signal = signal.Output_postprocess ;
subplot1(6) ; plot(tms, signal) ; axis('tight') ; title('$z(t)$ of KP-DS-MVDR','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;

filter = 'MVDR_DS_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output ; % signal = signal.Output_postprocess ;
subplot1(7) ; plot(tms, signal) ; axis('tight') ; title('$z(t)$ of KP-MVDR-DS','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; xticklabels({'0','1','2','3'}) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;
xlabel('time (s)') ;

a=findobj(gcf); % get the handles associated with the current figure

box on;

a=findobj(gcf); % get the handles associated with the current figure
allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Times New Roman','FontWeight','Normal','LineWidth',0.5,'FontSize',16);
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

linkaxes( allaxes , 'x' ) ;




figure(2); %1
subplot1( 7 , 1 , 'Min' , [0.05 0.07] , 'Max' , [0.8 1.01] , 'Gap' , [ 0.00 , 0.04] , 'XTickL' , 'Margin' , 'YTickL' , 'All' ) ; 

subplot1(1) ; plot(tms, SOI) ; axis('tight') ; title('$x(t)$','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;

subplot1(2) ; plot(tms, Input) ; axis('tight') ; title('$y_1(t)$','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;
       
filter = 'DS_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output_postprocess ;
subplot1(3) ; plot(tms, signal) ; axis('tight') ; title('$\hat z(t)$ of DS','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;

filter = 'MVDR_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output_postprocess ;
subplot1(4) ; plot(tms, signal) ; axis('tight') ; title('$\hat z(t)$ of MVDR','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;
ylabel('Amplitude') ;

filter = 'MVDR_K_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output_postprocess ;
subplot1(5) ; plot(tms, signal) ; axis('tight') ; title('$\hat z(t)$ of KP-MVDR','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;

filter = 'DS_MVDR_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output_postprocess ;
subplot1(6) ; plot(tms, signal) ; axis('tight') ; title('$\hat z(t)$ of KP-DS-MVDR','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;

filter = 'MVDR_DS_F' ;
signal = load([filters_dir,'/', filter , postname]) ; 
signal = signal.Output_postprocess ;
subplot1(7) ; plot(tms, signal) ; axis('tight') ; title('$\hat z(t)$ of KP-MVDR-DS','Position',[3.5,-0.5,1]) ;
xticks([0,1,2,3]) ; xticklabels({'0','1','2','3'}) ; yticks([-0.5,0.5]) ; yticklabels({'-0.5','0.5'}) ;
xlabel('time (s)') ;

a=findobj(gcf); % get the handles associated with the current figure

box on;

a=findobj(gcf); % get the handles associated with the current figure
allaxes=findall(a,'Type','axes');
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');

set(allaxes,'FontName','Times New Roman','FontWeight','Normal','LineWidth',0.5,'FontSize',16);
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

linkaxes( allaxes , 'x' ) ;