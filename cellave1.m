function [Y,N,N2,nr,nc,binlat,binlon,num,wnum,rr,cc] = cellave1(X,lonlat,latlim,lonlim,scale)
% cellave:  average site or station data into lat-lon cells
% CALL:  [Y,N,N2] = cellave1(X,lonlat,latlim,lonlim,scale);
%
%*********** IN **************
%
% X (mX x nX)r  station or site data, mX years, nX stations
% lonlat (nX x 2)r  lon, lat, in mapping decimal degree units
% latlim (1 x 2)r  latitude limits for the colormap
% lonlim (1 x 2)r   longitude limits for the color map
% scale (1 x 1)i  number of cells per degree in the color map
%
%*********** OUT *************
%
% Y (mY x nY)r  Time series matrix of cell-average values of X
% N (mY x 1)i  number of years with valid data in cell
% N2 (mY x nY)i number of sites in cell ave for each year
% [binlat binlon]  lats and lons of centers of cells with any entries
% [num wnum] number and areal-weighted number of 
%  rr,cc  row and col (s to n, w to e) of each site or station
%
%************ NOTES **********
%
% All time series in X must have data for all years.  In other
% words, no NaN values in X
%
% STEPS
%
% compute r,c and maplegend using sizem(latlim,lonlim,scale)
% determine row and col of each site using 
% make a dummy map using rand(r,c)
%		[rr,cc]=setpostn(map,maplegend,lonlat(:,2),lonlat(:,1));
% compute cell centers and number of sites in each cell using
%		[binlat,binlon,num,wnum]=histr(lonlat(:,2),lonlat(:,1));
% compute pointer matrix of cells to sites, and count # sites per cell
% Compute cell averages; assign NaN if no sites in cell


% compute number of rows and cols in lon-lat grid
[nr,nc,maplegend]=sizem(latlim,lonlim,scale);

% compute number of sites and years
[nyr,nsites]=size(X);

% Compute Logical pointer 3-d matrix from sites to cells
L=logical(zeros(nr,nc,nsites));

% Make a dummy map;
dummap=rand(nr,nc);
lats=lonlat(:,2); lons=lonlat(:,1);

% Compute row and col index for each site.  Rows go south to north,
% cols west to east
[rr,cc]=setpostn(dummap,maplegend,lats,lons);

% compute cell centers and number of sites for all cells having sites
[binlat,binlon,num,wnum]=histr(lats,lons,scale);
ncell1 = length(num);  % number of cells with any sites

% Make a logical index matrix tying sites to elements of L
Lcol = L(:);  % string L out as col vectornprod = nr*nc;
isite= (1:nsites)';
ind = sub2ind([nr,nc,nsites],rr,cc,isite); % index of site/row/col into Lcol
Lcol(ind)=logical(1);
L = reshape(Lcol,nr,nc,nsites);

% Compute sum of sites in each box
sum1 = sum(L,3); % 2-dim matrix of sums for each cell

%Initialize 
Y = repmat(NaN,[nr nc nyr]); % cell average
N = repmat(NaN,[nr nc]); % number of years of non-NaN data in cell
N2 = repmat(NaN,[nr nc nyr]); % number of sites going into cell ave in each year

% Compute cell means
for n1 = 1:nr;
   n1vect = repmat(n1,nsites,1);
   for n2 = 1:nc;
      indcol = sub2ind([nr nc],n1,n2); % linear index to cell row/col
      % will run  S to N along westmost, then next westmost, and so forth
      n2vect=repmat(n2,nsites,1);
      sum2 = sum1(n1,n2);
      if sum2==0;
         %Y(:,indcol)=NaN;
         %N(indcol)=NaN; % number of years of data (non-NaN) at the cell
         Y(n1,n2,:)=NaN;
         N(n1,n2)=NaN;
      else
         Lyes = rr==n1vect & cc==n2vect; % col vector pointing to sites
         X1 = X(:,Lyes); % grab the sub-tsm of data
         % if only one site, that is the cell average
         %N2(:,indcol)=0;
         N2(n1,n2,:)=0;
         if sum(Lyes)==1;
            %N2(~isnan(X1),indcol)=1;
            %Y(:,indcol)=X1;
            %N(indcol) =sum(~isnan(X1));
            N2(n1,n2,~isnan(X1))=1;
            Y(n1,n2,:)=X1;
            N(n1,n2)=sum(~isnan(X1));
            
         else
            X2=X1';
            %N2(:,indcol) = (sum(~isnan(X2)))';
            N2(n1,n2,:)=(sum(~isnan(X2)))';
            x2=nanmean(X2);
            Lvalid = ~isnan(x2);
            %Y(:,indcol)=x2';
            %N(indcol)=sum(Lvalid);
            Y(n1,n2,:)=x2';
            N(n1,n2)=sum(Lvalid);
         end
      end
   end
end
        
