S=zeros(2,length(x));
tpnt=zeros(1,4);

for j=1:1;       %j=1:length(i)
	p=0.0001;
	% detrend mean

	k1=1;  %  stay in the menu
	while k1==1
		k2=menu('SELECT ONE','INITIAL SPLINE FITS','QUIT')
		clc
		if k2==1;  % initial spline fits
			disp('WORKING ...')
			p1=1e-5; p2=1e-4;
			S(1,:)=csaps(t,x,p1,t);
			S(2,:)=csaps(t,x,p2,t);
			plot(t,x,t,S(1,:),t,S(2,:))
			v=axis;
			tpnt=[.1 * (v(2)-v(1))+v(1)   .9*(v(4)-v(3))+v(3)...
			.1 * (v(2)-v(1))+v(1)   .85*(v(4)-v(3))+v(3)];
			text(tpnt(1),tpnt(2),['p = ',num2str(p1)]);
			text(tpnt(3),tpnt(4),['p = ',num2str(p2)]);
			title(['SPLINE TO RINGWIDTHS: ',N(j,:)])

			xlabel('YEAR')
			ylabel(' mm x 100')
	
			pause

			disp('first time')
		elseif k2==2
			disp('quitting')
			k1=2;
		end
	end
	disp('working ...')



	s=csaps(t,x,p,t);
	plot(t,x,t,s)
end
