function treeplt1(x,yr,axlims,YR)
%
% D Meko 3-27-96
%
% Plot a long tree-ring or other time series on multiple plots, with
% subset of years on each plot.  Designed for portrait orientation,
% two sets of axes on each page, line plots with asterisks at data points.
%
%
%************** IN ARGS *******************
%
% x (mx x 1)r time series
% yr (mx x 1)i years for x
% axlims(? x 2)i x axis limits (years) for each plot
% YR (? x 2)i years of data to put on plots for axlims


[mYR,nYR]=size(YR);


x=100*x; % convert tree-ring index to pct normal growth

close all
orient landscape

x1=min(x);
x2=max(x);

% Compute number of plot pages
npages = ceil(mYR/2);

for i = 1:npages
	jgo=i*2-1; % starting row of YR for this page
	jstop=jgo+1; % ending...
	figure(i);
	subplot(2,1,1)
	L1 = yr>=YR(jgo,1) & yr<=YR(jgo,2);
	xx=x(L1);
	t = yr(L1);
	 
	plot(t,xx,'*')
	axis([axlims(jgo,1) axlims(jgo,2)  x1 x2]);
	line(t,xx);
	title('Tree Growth, San Pedro Martir White Fir (Site PMO)')
	ylabel('% of Normal Growth')
	grid

	if ~(rem(mYR,2)==1 & i==npages)

		subplot(2,1,2)
		L1 = yr>=YR(jstop,1) & yr<=YR(jstop,2);
		xx=x(L1);
		t = yr(L1);
	
		
		plot(t,xx,'*')
		axis([axlims(jstop,1) axlims(jstop,2)  x1 x2]);
		line(t,xx);
		grid
	end
end
