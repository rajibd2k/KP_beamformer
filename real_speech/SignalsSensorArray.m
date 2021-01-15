% % Signals : SOI and Interferences at the first sensor
% % theta : angle (in degrees) of the source signal and interferences impinging on the ULA
% % FS : sampling frequency of the signals
% % M : number of sensors
% % delta : inter sensor distance
% % c : velocity of the signals in the medium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Signals_ULA = SignalsSensorArray(Signals, theta, FS, M, delta, c)

num_signals = size(Signals,2) ;
num_samples = size(Signals,1) ;

theta = theta/180*pi ; 

Signals_ULA = zeros( num_samples , num_signals , M ) ;
Signals_ULA(:,:,1) = Signals ;

Ts = 1 / FS ;
t = [0 : (num_samples-1)]'*Ts ;

for idx_signal = 1 : num_signals
    
    Signal = Signals( : , idx_signal ) ;
%     gradient = 2 ./ (max(Signal) - min(Signal)) ;
%     intercept = - (max(Signal) + min(Signal)) ./ (max(Signal) - min(Signal)) ;
%     Signal = gradient*Signal + intercept ;
    Signal = Signal / max(abs(Signal)) ;
    Signals_ULA(:,idx_signal,1) = Signal ; 
    
    Signal = timeseries( Signal , t ) ;
    theta_idx = theta( idx_signal ) ;
   
    for idx_M = [2 : M]
    
        delay_at_sensor = (idx_M - 1) * delta * cos(theta_idx) / c ;
        t_reqd = t - delay_at_sensor ;
        signal = resample( Signal , t_reqd ) ;
        signal = signal.Data ;
        signal( find(isnan(signal)) ) = 0 ;
        
%         gradient = 2 ./ (max(signal) - min(signal)) ;
%         intercept = - (max(signal) + min(signal)) ./ (max(signal) - min(signal)) ;
%         signal = gradient*signal + intercept ;

        signal = signal / max(abs(signal)) ;

        Signals_ULA(:,idx_signal,idx_M) = signal ; 
        
    end 
    
end

end


% function Signals_ULA = SignalsSensorArray(Signals, theta, FS, M, delta, c)
% 
% num_signals = size(Signals,2) ;
% num_samples = size(Signals,1) ;
% 
% theta = theta/180*pi ; 
% 
% Signals_ULA = zeros( num_samples , num_signals , M ) ;
% Signals_ULA(:,:,1) = Signals ;
% 
% Ts = 1 / FS ;
% t = [0 : (num_samples-1)]*Ts ;
% 
% for idx_signal = 1 : num_signals
%     
%     Signal = Signals( : , idx_signal ) ;
%     theta_idx = theta( idx_signal ) ;
%     smooth_original = abs([zeros(1000,1) ; Signal ; zeros(1000,1)]) ;
%     for idx_smooth = 1:10
%         smooth_original = smooth( smooth_original ,101 ) ;
%     end
%     [y_maxima , x_maxima , y_minima , x_minima ] = extrema(smooth_original) ;
%     [~ , y_maxima_indices] = sort(y_maxima, 'descend') ; 
%     x_maxima = x_maxima(y_maxima_indices) ; x_maxima = x_maxima(1:3) ; % top 3 maxima
%     [~ , y_minima_indices] = sort(y_minima, 'descend') ; % not ascend
%     x_minima = x_minima(y_minima_indices) ; x_minima = x_minima(1:3) ; % bottom 3 minima
%     indices_original = [x_maxima ; x_minima] ;
%     
%     for idx_M = [2 : M]
%     
%         delay_at_sensor = (idx_M - 1)*delta * cos(theta_idx) / c ;
%         delay_type = sign(delay_at_sensor) ; % -1 : ahead , 1 : lag
%         
%         delay_samples = round( abs(delay_at_sensor) / Ts ) ; % samples shift at FS 
%         Ts_reqd = rem( abs(delay_at_sensor) , Ts ) ;
% 
%         if abs(delay_at_sensor) == 0
%             Signals_ULA(:,idx_signal,idx_M) = Signal ;
%             
% 
%         elseif and( rem( abs(delay_at_sensor) , Ts ) == 0 , delay_type == 1 )
%             sample_shift = round( abs(delay_at_sensor) / Ts ) ;
%             signal = [ zeros(sample_shift,1) ; Signal(1:end-sample_shift) ] ;
%             Signals_ULA(:,idx_signal,idx_M) = signal ;
%             
%         elseif and( rem( abs(delay_at_sensor) , Ts ) == 0 , delay_type == -1 )
%             sample_shift = round( abs(delay_at_sensor) / Ts ) ;
%             signal = [ Signal(sample_shift+1:end) ; zeros(sample_shift,1) ] ;
%             Signals_ULA(:,idx_signal,idx_M) = signal ;
%             
%         elseif and( Ts_reqd > 0 , delay_type == 1 )
%             p = round( 5 * Ts / Ts_reqd ) ; q = 1 ; % oversampling by 5 times than required
%             signal = resample(Signal,p,q) ; Ts_resampled = Ts / p ;
%             sample_shift = round( abs(delay_at_sensor) / Ts_resampled ) ;
%             signal = [ zeros(sample_shift,1) ; signal(1:end-sample_shift)] ;
%             signal = resample(signal,q,p) ;
%             
%             smooth_signal = abs([zeros(1000,1) ; signal ; zeros(1000,1)]) ;
%             for idx_smooth = 1:10
%                 smooth_signal = smooth( smooth_signal ,101 ) ;
%             end
%             [y_maxima , x_maxima , y_minima , x_minima ] = extrema(smooth_signal) ;
%             [~ , y_maxima_indices] = sort(y_maxima, 'descend') ; 
%             x_maxima = x_maxima(y_maxima_indices) ; x_maxima = x_maxima(1:3) ; % top 3 maxima
%             [~ , y_minima_indices] = sort(y_minima, 'descend') ; % not ascend
%             x_minima = x_minima(y_minima_indices) ; x_minima = x_minima(1:3) ; % bottom 3 minima
%             indices_signal = [x_maxima ; x_minima] ;
%             
%             diff_index = round( median(indices_signal - indices_original) - delay_samples ) ;
%             if diff_index >= 0
%                 signal = [ signal(diff_index + 1 : end) ; zeros(diff_index,1) ] ;
%             else
%                 signal = [ zeros( abs(diff_index) , 1) ; signal( 1 : end - abs(diff_index) ) ] ;
%             end
%             Signals_ULA(:,idx_signal,idx_M) = signal ;
%             
%         elseif and( Ts_reqd > 0 , delay_type == -1 )
%             p = round( 5 * Ts / Ts_reqd ) ; q = 1 ; % oversampling by 5 times than required
%             signal = resample(Signal,p,q) ; Ts_resampled = Ts / p ;
%             sample_shift = round( abs(delay_at_sensor) / Ts_resampled ) ;
%             signal = [ zeros(sample_shift,1) ; signal(1:end-sample_shift)] ;
%             signal = resample(signal,q,p) ;
%             
%             smooth_signal = abs([zeros(1000,1) ; signal ; zeros(1000,1)]) ;
%             for idx_smooth = 1:10
%                 smooth_signal = smooth( smooth_signal ,101 ) ;
%             end
%             [y_maxima , x_maxima , y_minima , x_minima ] = extrema(smooth_signal) ;
%             [~ , y_maxima_indices] = sort(y_maxima, 'descend') ; 
%             x_maxima = x_maxima(y_maxima_indices) ; x_maxima = x_maxima(1:3) ; % top 3 maxima
%             [~ , y_minima_indices] = sort(y_minima, 'descend') ; % not ascend
%             x_minima = x_minima(y_minima_indices) ; x_minima = x_minima(1:3) ; % bottom 3 minima
%             indices_signal = [x_maxima ; x_minima] ;
%             
%             diff_index = round( median(indices_original - indices_signal) - delay_samples ) ;
%             if diff_index >= 0
%                 signal = [ zeros(diff_index,1) ; signal(1 : end-diff_index) ] ;
%             else
%                 signal = [ signal( abs(diff_index) + 1 : end) ; zeros( abs(diff_index) , 1) ] ;
%             end
%             Signals_ULA(:,idx_signal,idx_M) = signal ;
%             
%         end
%        
%     end 
%     
% end
% 
% 
% for idx_M = 1 : M 
%     
%     Signals_max = ones(num_samples,1) * max(Signals_ULA(:,:,idx_M)) ;
%     Signals_min = ones(num_samples,1) * min(Signals_ULA(:,:,idx_M)) ;
%     m = 2 ./ (Signals_max - Signals_min) ;
%     c = - (Signals_max + Signals_min) ./ (Signals_max - Signals_min) ;
%     
%     Signals_ULA(:,:,idx_M) = m.*Signals_ULA(:,:,idx_M) + c ; 
%     
% end
% 
% end
