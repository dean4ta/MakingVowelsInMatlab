% Whisper or Sound order of:  [original AH EH EE OH OO]
function total = OneForWhisperTwoForSound(OneForWhisperTwoForSound, keynum)
load('AH A and B.mat')
load('EH A and B.mat')
load('EE A and B.mat')
load('OHa A and B.mat')
load('OOa A and B.mat')

dur = 1;
fs = 44100;
tt = 0:1/fs:dur;

if OneForWhisperTwoForSound == 1
    xx = rand(1,length(tt));
    %plot(xx)
elseif OneForWhisperTwoForSound == 2
    A = 1;
    f0 = 220*(2^((keynum-49)/12));
    xx = 0;
    for i = 1:30
        xx = xx + A*cos(i*f0*2*pi*tt);
    end    
else
    disp('please enter 1 Wisper or 2 for Sound')
end
% plot(xx)
%making sure the amplitudes are the same
xx = xx./max(xx);
AH = filter(B_ah, A_ah, xx);
AH = AH./max(AH);
EH = filter(B_eh, A_eh, xx);
EH = EH./max(EH);
EE = filter(B_ee, A_ee, xx);
EE = EE./max(EE);
OH = filter(B_oo, A_oo, xx);
OH = OH./max(OH);
OO = filter(B_oh, A_oh, xx);
OO = OO./max(OO);
silence = zeros(1,8000);
silencex = [silence silence silence silence silence];

total = [xx silencex AH silence EH silence EE silence OH silence OO];
soundsc(total, fs);


end