 function [outsig,noise,dbnoise]=addnoise(insig,orgnoise,rnl)

%disp('------------------ ADD STATIONARY NOISE ------------------');

insignal=insig';

%generate noise sequence using randn func

%orgnoise=rand(length(insig'),1);
	
%Input signal energy

insigenr=sum(insig'.*insig')/(length(insig'));

% size(insigenr)
%Orginal noise energy

orgnoiseenr=sum(orgnoise.*orgnoise)/(length(orgnoise));

%Required noise energy 

%    SNR=10 LOG[SP/NP]   => SP/NP=1O^(SNR/10)  LET IT BE  X
%    NP=SP/X 

outsig=zeros(length(rnl),length(insig));
noise=zeros(length(rnl),length(insig));

for i=1:length(rnl)
    
    y(i)=10^(rnl(i)/10);
    
    requirednoiseenr(i)=insigenr/y(i);

    %multiplication factor to be added for noise

    multfact(i)=sqrt(requirednoiseenr(i)/orgnoiseenr);

    %final noise to be mutliplied

    noise(i,:)=orgnoise*multfact(i);
    dbnoise(i)=10*log10(sum(noise(i,:)'.*noise(i,:)')/(length(noise(i,:)')));

    %output signal with noise

    outsig(i,:)=insig'+noise(i,:);
    outsig(i,:)=outsig(i,:)./max(abs(outsig(i,:)));
       
    %final noise energy
    finalnoiseenr(i)=sum(noise(i,:).*noise(i,:))/(length(noise(i,:)));

    %final SNR value
    Addedsnr(i)=10*log10(insigenr/finalnoiseenr(i));
end

