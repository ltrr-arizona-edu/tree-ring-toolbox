function daychk1(day,yr,xp,xe,xeold)
% View and check daily ppt estimation yielded by dayest.m
%
% D Meko   11-9-93
%
%
%**********   PRELIMS ***********
%
% In main workspace, load xdaily.dat, and assign appropriate columns
% as   yr, day, xk, xp  (see dayest2.m)
% In main workspace, set  xeold=xk
% In main workspace, edit dayest2.m to restrict julian-day range, and
%   run to get xe
%
%



k=1;


while k==1;
d = input('Julian Day to check out: ');
L=day==d;
plot(yr(L),xeold(L),yr(L),xe(L),'*',yr(L),xp(L),'+');
title(['Day ',int2str(d),' ,   X is estimate']);
[yr(L)  xp(L)  xeold(L)  xe(L)]
k=input('Another day? [1/0] ');
end

