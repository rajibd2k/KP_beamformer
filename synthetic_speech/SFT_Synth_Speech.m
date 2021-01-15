% Generatic Synthetic Speech using Formant Synthesizer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dur : duration of the speech signal to be generated in seconds
% FS : sampling frequency in Hz
% F0 : pitch frequency in Hz
% F : formant frequency vector in Hz , F = [F1;F2;F3;F4;F5]
% Bw : formant bandwidth vector in Hz , Bw = [Bw1;Bw2;Bw3;Bw4;Bw5]
% speech_type : 'voiced' / 'unvoiced'

function [ excitation , glottal , v1 , v2 , v3 , v4 , v5 , vts , lips , speech ] = SFT_Synth_Speech( dur , FS , F0 , F , Bw , speech_type  )

nsamps = floor(dur*FS) ;     % number of samples

%***********************************************************************************************************************
% Vocal Tract System
%***********************************************************************************************************************

% input analog frequencies corresponding to the poles
F1 = F(1) ; F2 = F(2) ; F3 = F(3) ; F4 = F(4) ; F5 = F(5) ;

% input analog bandwidths corresponding to the poles
Bw1 = Bw(1) ; Bw2 = Bw(2) ; Bw3 = Bw(3) ; Bw4 = Bw(4) ; Bw5 = Bw(5) ;

% number of actual formants
num_F = max( find( F > 0 ) ) ;

% pole radii
R1 = 1 - Bw1/FS ; R2 = 1 - Bw2/FS ; R3 = 1 - Bw3/FS ; R4 = 1 - Bw4/FS ; R5 = 1 - Bw5/FS ;

% digital frequencies corresponding to the poles
w1 = 2*pi*F1/FS ; 
w2 = 2*pi*F2/FS ; 
w3 = 2*pi*F3/FS ; 
w4 = 2*pi*F4/FS ; 
w5 = 2*pi*F5/FS ; 

% poles
if strcmp( speech_type , 'voiced')
    poles_glottal = [ 0.98 ; 0.98 ] ; % voiced
else
    poles_glottal = [ 0 ] ; % unvoiced
end

poles_v1 = R1 .* exp(j*w1) ; poles_v1 = [ poles_v1 ; conj( poles_v1 ) ] ; 
poles_v2 = R2 .* exp(j*w2) ; poles_v2 = [ poles_v2 ; conj( poles_v2 ) ] ; 
poles_v3 = R3 .* exp(j*w3) ; poles_v3 = [ poles_v3 ; conj( poles_v3 ) ] ; 
poles_v4 = R4 .* exp(j*w4) ; poles_v4 = [ poles_v4 ; conj( poles_v4 ) ] ;
poles_v5 = R5 .* exp(j*w5) ; poles_v5 = [ poles_v5 ; conj( poles_v5 ) ] ;
poles_vts = [poles_v1 ; poles_v2 ; poles_v3 ; poles_v4 ; poles_v5] ;    
poles_vts = poles_vts( 1 : 2*num_F ) ; % select the number of resonant poles

poles_lips = 0 ;

poles = [ poles_glottal ; poles_vts ; poles_lips ] ; 

% zeros
zer0s_glottal = zeros( size(poles_glottal) ) ;

zer0s_v1 = zeros( size(poles_v1) ) ;
zer0s_v2 = zeros( size(poles_v2) ) ;
zer0s_v3 = zeros( size(poles_v3) ) ;
zer0s_v4 = zeros( size(poles_v4) ) ;
zer0s_v5 = zeros( size(poles_v5) ) ;
zer0s_vts = zeros( size(poles_vts) ) ;

zer0s_lips = [ 0.98 ] ;

zer0s = [ zer0s_glottal ; zer0s_vts ; zer0s_lips ] ; 

% Transfer functions
[B_glottal,A_glottal] = zp2tf( zer0s_glottal , poles_glottal , 1) ; % control

[B_v1,A_v1] = zp2tf( zer0s_v1 , poles_v1 , 1) ; % control
[B_v2,A_v2] = zp2tf( zer0s_v2 , poles_v2 , 1) ; % control
[B_v3,A_v3] = zp2tf( zer0s_v3 , poles_v3 , 1) ; % control
[B_v4,A_v4] = zp2tf( zer0s_v4 , poles_v4 , 1) ; % control
[B_v5,A_v5] = zp2tf( zer0s_v5 , poles_v5 , 1) ; % control
[B_vts,A_vts] = zp2tf( zer0s_vts , poles_vts , 1) ; % control

[B_lips,A_lips] = zp2tf( zer0s_lips , poles_lips , 1) ; % control

[B,A] = zp2tf( zer0s , poles , 1) ; % control


% Impulse Train excitation
%------------------------------------------------------   
periodicity = round(  FS/F0 ) ;   % FS/F0 must be an integer , equal to the time period of the impulse train
F0_refined = FS/periodicity ;

if strcmp( speech_type , 'voiced')
    excitation = zeros( nsamps , 1 ) ;
    excitation( [1 : periodicity : nsamps]' )  = 1 ;
else
    excitation =  random( 'Normal' , 0 , 1 , [nsamps , 1] ) ;
end

% Generation of Synthetic Speech and its components
%------------------------------------------------------
glottal = filter( B_glottal, A_glottal ,excitation) ; %glottal = glottal - mean( glottal ) ; glottal = glottal / max( abs( glottal ) ) ;

v1 = filter( B_v1, A_v1 ,excitation) ; %v1 = v1 - mean( v1 ) ; v1 = v1 / max( abs( v1 ) ) ;
v2 = filter( B_v2, A_v2 ,excitation) ; %v2 = v2 - mean( v2 ) ; v2 = v2 / max( abs( v2 ) ) ;
v3 = filter( B_v3, A_v3 ,excitation) ; %v3 = v3 - mean( v3 ) ; v3 = v3 / max( abs( v3 ) ) ;
v4 = filter( B_v4, A_v4 ,excitation) ; %v4 = v4 - mean( v4 ) ; v4 = v4 / max( abs( v4 ) ) ;
v5 = filter( B_v5, A_v5 ,excitation) ; %v5 = v5 - mean( v5 ) ; v5 = v5 / max( abs( v5 ) ) ;
vts = filter( B_vts, A_vts ,excitation) ; %vts = vts - mean( vts ) ; vts = vts / max( abs( vts ) ) ;

lips = filter( B_lips, A_lips ,excitation) ; %lips = lips - mean( lips ) ; lips = lips / max( abs( lips ) ) ;

speech = filter( B, A ,excitation) ; %speech = speech - mean( speech ) ; speech = speech / max( abs( speech ) ) ;

