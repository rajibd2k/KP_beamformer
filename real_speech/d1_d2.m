%%% Find virtual and actual steering vectors 
%%% Delta : inter-sensor positions of the original ULA of M = M1 x M2 sensors
%%% Delta1, Delta2 : inter-sensor positions of the virtual arrays
%%% theta : Angle of incidence of the signal in degrees
%%% f : frequency range betweeen [0 1)
%%% c : speed of plane wave propagation (320 m/s)
%%% Ts : sampling period of the signals (32 kHz)
%%% d1, d2, d : Steering vectors of the signal at each frequency, using the virtual arrays

function [ Steering_vectors , Sensor_positions ] = d1_d2( Delta, M2, theta, f, c, Ts )

    M = length( Delta ) ;
    theta = theta/180*pi ;
    
    Delta2 = Delta(1:M2) ;
    
    M1 = M / M2 ;
    delta1 = Delta(M2+1) ;
    Delta1 = delta1*[0:(M1-1)]';
    
    Sensor_positions.Delta = Delta ;
    Sensor_positions.Delta1 = Delta1 ;
    Sensor_positions.Delta2 = Delta2 ;
    
    % steering vector of the ULA
    [f_mat , m_mat] = meshgrid(f, Delta/(c*Ts) ) ;
    d = exp(-1i*2*pi*f_mat.*m_mat*cos(theta)) ;
    
    % steering vector of the second virtual ULA
    [f_mat , m_mat] = meshgrid(f, Delta2/(c*Ts) ) ;
    d_2 = exp(-1i*2*pi*f_mat.*m_mat*cos(theta)) ;
    
    % steering vector of the first virtual ULA
    [f_mat , m_mat] = meshgrid(f, Delta1/(c*Ts) ) ;
    d_1 = exp(-1i*2*pi*f_mat.*m_mat*cos(theta)) ;
    
    Steering_vectors.d = d ;
    Steering_vectors.d_1 = d_1 ;
    Steering_vectors.d_2 = d_2 ;

end


