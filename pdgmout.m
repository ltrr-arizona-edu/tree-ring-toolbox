% pdgmout.m 

% Runs on output from pdgm1.m.  For getting output arrays for plotting
% in Splus


windx1=ftemp+xover;
windx2=fline+xover;
windy1=wsave .* yfact + yover;
windy2=straight+yover;


www1=[ff' ysm1 ysm2 ysmlo ysmhi];
www2=[windx1' windy1];
www3=[windx2' windy2];

save pdgm1.dat www1 -ascii
save pdgm2.dat www2 -ascii
save pdgm3.dat www3 -ascii
