function [rbt, D, W]=stats2(IT,Tnms,ITyrs)
% stats2: statistics on tree indices
% [rbt, D, W]=stats2(IT,Tnms,ITyrs);
% Last revised 9-2-99
%
% In tree-ring standardization, computes statistics on tree indices
%
%*** INPUT
%
% IT (? x 1) - Tree indices 
% Tnms (ns x 6) - Tree names
% Tyrs (ns x 3) - Year matrix
%
% OUTPUTS : 
%	rbt (1 x 2) - Average between-tree correl, and number
%		of correlations it is based on.
%	W (ns x 2) Wigley statistics -- EPS in col 1, SSS in col 2
%	D (ns x 6) basic statistics on tree indices; by column:
%		1 sample size (number of Non-NaN years)
%		2 mean
%		3 std dev
%		4 mean sensitivity
%		5 first-order autocorrel coef
%		6 95% conf lim for first-order autocorr coef
%		7 second-order parial ac coef
%
%*** UW FUNCTIONS CALLED
%
% c:\mlb\pacf
% meansen1
% rtree3
% c:\mlb\acf
% basic1
% wigley1
% rtree4
%
%*** TOOLBOXES NEEDED
% System Identification
%
% DESCRIPTION : [avcor,eps,sss,ms,r1,phi]=stats2(IT,Tnms,Tyrs);
% Returns statistics on tree indices
%


% Number of lags for acf and pacf
nlags=3;

% Get the sample size
[NN,dum]=size(Tnms);


rr = zeros(NN,3);

% Find the correlation coefficients between trees
r=rtree3(IT,Tnms,ITyrs);

% Find the average correlation between trees
avcor=rtree4(r);

% Get the EPS and SSS statistics
[sss,epops]=wigley1(avcor(1),NN);

% Get sample size, mean and standard dev
[mnx,stdx,nyears] = basic1(ITyrs,IT);

% Get the mean sensitivity
ms = meansen1(IT,ITyrs);

% Get autocorrelation and partial autocorrelation coefficients
for scid=1:NN,
 % Cull out individual Tree indices
 yrv=(ITyrs(scid,1):ITyrs(scid,2))';
 xv=IT(ITyrs(scid,3):ITyrs(scid,3)+ITyrs(scid,2)-ITyrs(scid,1));

 % Check for segmentation
 lzn=isnan(xv);
 zn=xv;
 zn(lzn)=[];
 if ~isempty(zn),
  lzo=zeros(length(lzn),1);
  for k=1:length(lzn)-1,
    if lzn(k)~=lzn(k+1),
      lzo(k)=k;
    end
  end
  lz=lzo;
  lz(find(lzo==0))=[];

  if ~isnan(xv(1)),
    lz=[0;lz];
  end
  if ~isnan(xv(length(xv))),
    lz=[lz;length(xv)];
  end
   
  [r1,SE2,r95]=acf(zn,nlags);
  [phi,SE3] = pacf(zn,nlags);

 end
rr(scid,:) = [r1(1) r95 phi(2)];  % First-order acf and its 95 conf lim
end

rbt=avcor;
W = [epops sss] ; % Wigley et al info
D = [nyears mnx stdx ms rr];

% End of file

