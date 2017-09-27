% --------------play_fugue.m-------------- %
load bach_fugue.mat
fs = 44100; % 11025 Hz also works
%theVoices(1).durations = theVoices(1).durations;
sp = 0.125;
xx = zeros(1, ceil(sum(theVoices(1).noteNumbers)*sp*fs));

for i = 1:length(theVoices)
for kk = 1:length(theVoices(i).noteNumbers)
    keynum = theVoices(i).noteNumbers(kk);
    tone = glottalkey2note(keynum, theVoices(i).durations(kk)*sp);% <------- Fill in this line
    strt = theVoices(i).startPulses(kk)*fs*sp;
    lend = strt+length(tone)-1;
    xx(strt:lend) = xx(strt:lend) + tone;
end
end
%soundsc(x1, fs)
% 
% theVoices(2).durations = theVoices(2).durations;
% x2 = 0;
% for kk = 1:length(theVoices(2).noteNumbers)
%     keynum = theVoices(2).noteNumbers(kk);
%     tone = glottalkey2note(keynum, theVoices(2).durations(kk));% <------- Fill in this line
%     x2 = [x2 tone];
% end
% %soundsc(x2, fs)
% 
% theVoices(3).durations = theVoices(3).durations;
% x3 = 0;
% for kk = 1:length(theVoices(3).noteNumbers)
%     keynum = theVoices(3).noteNumbers(kk);
%     tone = glottalkey2note(keynum, theVoices(3).durations(kk));% <------- Fill in this line
%     x3 = [x3 tone];
% end
% %soundsc(x3, fs)
% 
% x2 = [zeros(1, length(x1)-length(x2)), x2];
% x3 = [zeros(1, length(x1)-length(x3)), x3];
% xx = x1+x2+x3;
soundsc(xx, fs)