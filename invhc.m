
%function invhc
% sets default inverthardcopy for figures on or off
% Will want inverthardcopy off if 
%	Plots have colored or white lines on black bkgd, and
%	Want color slides

txt1='defaultfigureinverthardcopy';


setoff = ['set(0,''defaultfigureinverthardcopy'',''off'')'];
seton  = ['set(0,''defaultfigureinverthardcopy'',''on'')'];


k = input('Set Default InvertHardCopy to off?  [Y]','s');


if isempty(k)  | k=='Y'  |k=='y'
	k='Y';
else
	k='N';
end


if k=='Y'
	eval(setoff)
else
	eval(seton)
end
