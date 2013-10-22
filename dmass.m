
% dmass.m   double mass plots of monthly total ppt series

%*****************   PRELIMS   ********************************

% Run sequence FIL999.F, rear1.m to prepare arrays in Y.MAT
% Prepare an array E, ?x10, specifying the dmass plots

%******************* PRELOADS  **************************

% Y.MAT
% E.DAT

%********************************************************

% Temporary initializng

E=[8 3 1 5 3   0  0  0  0  0;
10 4 1 5 11 3  0  0  0  0];

mtot=130;  % total number of years in axis


%*******************************************************

Z = zeros(mtot,4);
Y = zeros(mtot,4);
X = zeros (mtot,4);

k3=0;


while k3~=4

k3=menu('PICK ONE','NEW TARGET SERIES','SAME SERIES, DIFF MONTH',...
'REPEAT SAME SERIES AND MONTH','QUIT');


if k3==1, k1 = input('TARGET SERIES NUMBER:  '); , end
if k3==2 | k3==1,  k2 = input('MONTH [1,2,...] : '); , end
if k3==4, break, end

d = E(k1,:); % rv of target series number, num of pairs, and which

nump = sum(d~=0) - 2;


for k=1:nump  % loop for each pair
	disp('here')
	ser= d(k+2);
	j1=(d(1)-1)*12+k2;
	j2=(ser-1)*12+k2;
	T=[b   D(:,j1)   D(:,j2)];
	L=T(:,2:3)>=99.90;
	l=(any(L'))';
	T(l,:)=[];
	[m1,n1]=size(T);
	M1(k)=m1;

	T=[T(:,1)   cumsum(T(:,2))  cumsum(T(:,3))];

	X(:,k)= [T(:,1);zeros(mtot-m1,1)];
	Y(:,k)= [T(:,2);zeros(mtot-m1,1)];
	Z(:,k)=[T(:,3); zeros(mtot-m1,1)];
end

for k=1:nump;

	x=X(1:M1(k),k);
	y = Y(1:M1(k),k);
	z = Z(1:M1(k),k);

	xx=x;
	num = floor(M1(k)/10);
	xxx=xx(1:10:M1(k));
	yyy=y(1:10:M1(k));
	zzz=z(1:10:M1(k));

	nm=220+k;
	subplot(nm);
	plot(Z(1:M1(k),k),Y(1:M1(k),k),'+',zzz,yyy,'o');
	xlabel(F(d(k+2),1:6));
	ylabel(F(d(1),1:6));
	text(.4,1.0,['MONTH ',num2str(k2)],'sc')

	xshort=rem(xxx,100);
	for i=1:length(xxx)
		 text(zzz(i),yyy(i),int2str(xshort(i)));
	end
	pause
end

	k=[];
	k=input('TIME-SERIES PLOTS? Y/N [N]: ','s');
		if isempty(k), k='N'; end;
		if k=='Y'
			subplot 
			for k=1:nump
				sx= X(1:M1(k),k);
				dz = diff(Z(1:M1(k),k));
				dy = diff(Y(1:M1(k),k));
				sz = [Z(1,k); dz];
				sy = [Y(1,k); dy];
				plot(sx,sy,'-',sx,sz,'--')

				text(.2,.90,['SOLID -- ',F(d(1),:)],'sc')
				text(.2,.85,['DASHED -- ',F(d(k+2),:)],'sc')
				text(.7,.90,['MONTH -- ',num2str(k2)],'sc')
			
				r= corrcoef(sz,sy);
				r=r(1,2);
				text(.7,.85,['CORRELATION = ',num2str(r)],'sc')				

				xlabel ('YEAR')
				ylabel ('PPT (INCHES)')
				pause
		   end
		end
			k=[];
			k = input('ADJUST THE TARGET SERIES? Y/N [N]','s');
			if isempty(k), k='N'; end
			if k=='Y'
				k7=[];
				k7= input('PRINT THE YEARS TABLE? Y/N  [N]','s');
				if isempty(k7), k7='N';
				end
				if k7=='Y'
					disp('PRINT SCREEN OF STATION NAMES');
					pause(2)
					clc, home;
					disp([int2str(d(1)),' ',F(d(1),:)]);
					disp(' ')
					for k8 = 1:nump
						
						disp([int2str(k8),' -- ',...
						int2str(d(k8+2)),' ',F(d(k8+2),:)]);
					end   % of k8 loop
					disp(' ')
					disp('MONTH')
					disp(k2)
					disp(' ')
					disp('SHFT-PRTSCRN NOW')
					pause
					clc
					disp(' ')
					JJ=X~=0;
					JJ=sum(JJ);
					mx=max(JJ);
					XT=X(1:mx,:);
					save XT.dat XT -ascii
%					!print XT.dat
					 
				end   % end of k7 loop for printing the years table
		
				k7 = input('PLOT FRAME OF DIAGNOSTIC PAIR?  1-4:  ')
				subplot(111);
				k=k7;



				x=X(1:M1(k),k);
				y = Y(1:M1(k),k);
				z = Z(1:M1(k),k);

				xx=x;
				num = floor(M1(k)/10);
				xxx=xx(1:10:M1(k));
				yyy=y(1:10:M1(k));
				zzz=z(1:10:M1(k));

	         plot(Z(1:M1(k),k),Y(1:M1(k),k),'+',zzz,yyy,'o');
				xlabel(F(d(k+2),1:6));
				ylabel(F(d(1),1:6));
				text(.4,1.0,['MONTH ',num2str(k2)],'sc')

				xshort=rem(xxx,100);
				for i=1:length(xxx)
		 			text(zzz(i),yyy(i),int2str(xshort(i)));
				end

				pause

% ANOVA test on adjacent segments of key pair
% Invokes Meko function vrat.m, based on variance between vs within
%   samples.

				dy = diff(Y(1:M1(k),k));
				sy = [Y(1,k); dy];
				I1=find(x>=YRS(1,1) & x<=YRS(1,2)); 
				I2=find(x>=YRS(2,1) & x<=YRS(2,2));
				[vr,dfb,dfw] = vrat(sy(I1),sy(I2));

				disp('F-RATIO   DF-BETWEEN  DF-WITHIN')
				disp(' ')
				disp([vr dfb dfw])
				pause
				clc
				home 
				

% Decide whether to adjust record.

				k7=menu('CHOOSE ONE','ADJUST EARLY SEGMENT',...
				'ADJUST LATER SEGMENT','DONT ADJUST');

				if k7==1   % Adjust the early segment
					tz1=z(I2);
					ty1=y(I2);
					tz2=z(I1);
					ty2=y(I1);
				elseif k7==2;  % Adjust later segment
					tz1=z(I1);
					ty1=y(I1);
					tz2=z(I2);
					ty2=y(I2);
				elseif k7==3;   % Dont continue with adjustment.
				end
				[cc,dd,ratio]=slp1(tz1,ty1,tz2,ty2);
				newy=sy;
				if k7==2
					ytemp=sy(I2) * ratio;
					newy(I2) = ytemp;
				elseif k7==1
					ytemp=sy(I1) * ratio;
					newy(I1)= ytemp
				end


				dz = diff(Z(1:M1(k),k));
				sz = [Z(1,k); dz];
				plot(X(I2,k),ytemp,X(1:M1(k),k),sy,X(1:M1(k),k),sz)
				title('ADJUSTED AND ORIGINAL')
				pause

%  Note that in both cases, t1 will be adjusted 

		end     % of loop for adjusting series
clg
subplot

end  % of while loop
