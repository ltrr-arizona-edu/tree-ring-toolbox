% mvgm.m  plots sample and model variograms for exponential and
%         spherical models, and sample variogram data.
%       
% For equations, see p. 374 in Isaaks and Srivastava, 1989
%
% Before using this pgm, you must run GEOEAS and model the 
% variogram(s).  Then must key the data from geoeas into an array
% with col 1 a sequence, 2 a ?, 3 the distance, and 4,5,... the
% sample variogram values.  Also need respective nugget, sill etc

% Will need to tailor to suit whatever cols hold the data you
% want to analyze.

% Set up to plot vgms for two maps on same plot, with sample values
% and models.

%*******************  PRELOAD  ********************************

% vgms.dat - holding data as described above

%************************************************************** 


clg;
h=vgms(:,3);
g2=vgms(:,5);  % year 1958, example for exponential vgm model
g4=vgms(:,7);  % year 1977, example for spherical vgm model

c1=1.5;  %constant needed later
c2=0.5;  % constant

n2=0.009;  % nugget
s2=0.020;  %sill
r2=160;  % range
m=length(h);

nn2=n2(ones(m,1),:);
ss2=s2(ones(m,1),:);

y2=nn2+ s2 * (ones(m,1) - exp (-(3/r2) .* h));


n4=0.02;
s4=0.023;
r4=90;

nn4=n4(ones(m,1),:);
ss4=s4(ones(m,1),:);

b=h /r4;
bb = b .* b .* b;

y4=nn4 + s4 * ((c1(ones(m,1),:) / r4) .* h  -  ...
   c2(ones(m,1),:) .* bb);



% set model sph vgm equal to sill for h> range

I=h>r4;
y4(I)=s4(ones(sum(I),1),:) + n4;


axis([0 500 0 .08]);
plot(h,g2,'*r',h,y2,'-r',h,g4,'+g',h,y4,'-g');
xlabel ('DISTANCE (KM)');
ylabel ('GAMMA');
text(152.1,0.025,'1958');
text(119.5,0.044,'1977');

