function u=oedeconv(Z,yr,y,yry,nn)
% oedeconv:  deconvolve climate input from long-term tree-ring series fit by OE model
% CALL: u=oedeconv(Z,yr,y,yry,nn);
%
% Meko 8-18-97
%
%************ IN ***************
%
% Z (mZ x 2)r time series of output (col 1) and input (col 2) used to fit model
% yr (mZ x 1)i  year vector for Z
% y (my x 1)r  long-term tree-ring index (long version of col 1 of Z);
% yry (my x 1)i year vector for y
% nn (1 x 3)i   order of OE model used to fit Z. B,F,and k orders given by 
%		nn(1,2,3).


% Remove Z-period mean from Z and y
meany = mean(Z(:,1));
meanu = mean(Z(:,2));

Z1 = [Z(:,1)-meany  Z(:,2)-meanu];
y1 = y - meany;


% Set model aliases
if nn(1)==1 & nn(2)==0; % OE(1,0)
   modtype = 1;
elseif nn(1)==2 & nn(2)==0; % OE(2,0)
   modtype = 2;
elseif nn(1)==1 & nn(2)==1; % OE(1,1)
   modtype = 3;
else
end

% Estimate parameters of OE model
th = oe(Z1,nn);


% Get the model parameters
switch modtype;
case 1;
   B=th(3,1);
   D=B;
case 2;
   B = th(3,1:2);
   D = B;  % polynomial to divide tree-ring series by
case 3; % OE(1,1)
   
otherwise
   error('only handles OE(1,0) or OE(2,0)so far');
end


% Deconvolve the short climate signal from the tree-ring series
[Q,R]=deconv(Z1(:,1),D);
% Add back the mean
u2a = Q + meanu;
yru2a = yr(1:length(u2a));


% Deconvolve the long climate signal from the tree-ring series
[Q,R]=deconv(y1,D);
% Add back the mean
u2b = Q + meanu;
yru2b = yry(1:length(u2b));








u = 1;
