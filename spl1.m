function s=spl1(X,n,yrgo, yrstop)

% Fits spline curve with parameter p to n+1 th col of time series matrix X
% Col 1 of X is assumed to hold the years
% Curve is fit to segment specified by [yrgo yrstop]

L=X(:,1)>=yrgo & X(:,1)<=yrstop;  % logical pointer to selected years
t=X(L,1);  % year column
x=X(L,n);  % the core to fit


p=input('INITIAL SETTING FOR SLINE STIFFNESS PARAMETER P: ');
s=csaps(t,x,p,t);
plot(t,x,t,s,'-');
jreset=1;
pltext(.1,.9,['p = ',num2str(p)])
jreset=j+1;
hold on

pause

j=1;  % while j ne 0, keep adding curves

while j==1;
	k1=menu('Choose one','Erase splines','Try another spline',...
		'Save plot as meta file', 'Quit')

	if k1==1;  %Erase splines from plot and start over
		clg
		jset=1
		hold off
	   p=input('INITIAL SETTING FOR SLINE STIFFNESS PARAMETER P: ');
		s=csaps(t,x,p,t);
		plot(t,x,t,s,'-');
		pltext(0.1,1.0-jset*0.1,['p= ',num2str(p)])
		hold on
		pause

	elseif k1==2;  %  Try another spline

		p=input('WHAT VALUE FOR STIFFNESS PARAMETER P? ')
		s=csaps(t,x,p,t);
		plot(t,s,'g--')
		pltext(0.1,0.8,['p= ',num2str(p)])
		pause

	elseif k1==3;  % Save current plot as metafile

		keyboard

	elseif  k1==4;  % get out
		hold off
		clg
		j=0;
		break
	end;  % of if

end  ;  % of while


