% --------------play_lamb.m-------------- %
mary.keys = [44 42 40 42 44 44 44 42 42 42 44 47 47];
% NOTES: C D E F G
% Key #40 is middle-C
mary.durations = .25 * ones(1,length(mary.keys));
fs = 44100; % 11025 Hz also works
% xx = zeros(1, sum(mary.durations)*fs + length(mary.keys));
xx = glottalkey2note(mary.keys(1), mary.durations(1));
for kk = 2:length(mary.keys)
    keynum = mary.keys(kk);
    tone = glottalkey2note(keynum, mary.durations(kk));% <------- Fill in this line
    xx = [xx tone];
end
specgram(xx)
soundsc(xx, fs)
