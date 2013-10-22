function rwlook1

% Visual check of cross-dating of ring width for pair of cores
%
% D. Meko   12-8-93;   last revised 9-22-94
%
%
% ******************   USER FUNCTIONS NEEDED  ***************
%
% flistin.m
% quantile.m
% rwread.m
% rwchng.m
% kendtau.m
% sgntst.m
% binom1.m
%  
% signtbl.mat
% cona12.mat   -- Conover table A12 for testing tau statistic
%
%***********************************************************

% Load lookup table for signif test of Kendall's tau
if ~exist('cona12');
  load cona12;
end

fonty=12; % hard-coded font size for call to pltext.m

disp('Reading the ring width file');
[zname,zyrs,zvar,zindx,fl]=rwinp;
if fl==-1,   % Failure in reading a file breaks the loop
  break;
  return;
end

kmnu=1;
while kmnu~=2,
  kmnu=shlmnu('Choose one','Pick new Series ?','QUIT');
  if kmnu==1,
    [x,yx,Xfn]=rwpick(zname,zvar,zyrs,zindx,'X');
    X=[yx,x];
    [x,yx,Yfn]=rwpick(zname,zvar,zyrs,zindx,'Y');
    Y=[yx,x];

    % Set vectors with rw and years
    x=X(:,2);  % series
    yrx=X(:,1);  % year vector
    y=Y(:,2);  % series 2
    yry=Y(:,1);  % year vector
    nx=length(x);
    ny=length(y);

    % Transform ring width
    w=rwchng(x);
    z=rwchng(y);
    yrw=yrx(2:nx);
    yrz=yry(2:ny); 

    k1=1;        % while control for level-1 menu
    while k1~=4;  % if k1==4, will quit the analysis
      k1=menu('Choose one: ','tau',...
      'Zoom RW','Zoom D','Quit');

      if k1==1;   % Initial sign and kendal-tau tests and plots
        close;
        % do this only on transformed data
        len=input('Length of Segment: ');
        offset=input('Offset (yrs): ');
        % Find overlap between transformed rw series
        yrso=[max([yrw(1) yrz(1)]) min([yrw(nx-1) yrz(ny-1)])];
        yro=(yrso(1):yrso(length(yrso)))';
        L3w=yrw>=yrso(1) & yrw <=yrso(length(yrso));
        L3z=yrz>=yrso(1) & yrz <= yrso(length(yrso));

        [yreh,H]=pullseg1(w(L3w),yrw(L3w),len,offset); % for w
        [yreb,B]=pullseg1(z(L3z),yrz(L3z),len,offset); % for z
        [mH,nH]=size(H);
        v=zeros(nH,1);

        for i=1:nH;
          h=H(:,i);
          b=B(:,i);
          [T,tau]=kendtau(h,b);
          v(i)=tau;
        end

        % 0.95 and 0.99 one-sided confidence levels for tau
        % Revision 1-31-94
        % One-sided test appropriate.  H0 is that tau=0.  H1 is that
        %   tau>0.  Function not written to test for significance of
        %   negative correlation, because unlikely to be of interest in
        %   comparing ring-width series.

        Tcrit=table1(cona12,len); % get a row from the lookup table
        tau95=Tcrit(2)/(len*(len-1)/2);  % for 95% signif level;
        % Denominator converts T to tau
        tau99=Tcrit(4)/(len*(len-1)/2);  % 99% signif level
        figure('Color','k');
        plot(yreh,v,yreh,tau95(ones(nH,1),:),'m--',...
        yreh,tau99(ones(nH,1),:),'m--');
        title(['KENDALLS TAU:        ',Xfn,'  and  ',Yfn]);
        pltext(.1,.9,fonty,[num2str(len),'-year segments']);
        pltext(.1,.8,fonty, [' Offset by ',num2str(offset),' years']);
        ylabel('Kendalls tau');
        xlabel('Ending year');
        hh=uicontrol('Style','Pushbutton',...
          'Position',[.9 .9 .1 .1],'Units','Normalized',...
          'Callback','print -dps','String','Laser');

      elseif k1==2;  % Initial and Zoom plots on ring width
        figure('Color','k');
        plot(yrx,x,'-',yry,y,'--');
        title(['RW: Solid= ',Xfn,'    Dashed= ',Yfn]);
        xlabel('Year');
        ylabel('RW (mm)');

        k=1;
        while k~=4;
          [t,v]=ginput(2);  % get coordinates of end points of segment
          tbeg=ceil(t(1));
          tend=floor(t(2));
          if tbeg< max([yrx(1) yry(1)])
            tbeg=max([yrx(1) yry(1)]);
          end
          if tend > min([yrx(length(x))  yry(length(y))]);
            tend=min([yrx(length(x))  yry(length(y))]);
          end

          tseg=(tbeg:tend)';  % time segment for zoom
          ix = tseg-yrx(1)+1;  % index into x
          iy = tseg-yry(1)+1;  % index into y
          figure('Color','k');
          plot(tseg,x(ix),'-',tseg,y(iy),'--');
          title(['RW: Solid= ',Xfn,'    Dashed = ',Yfn]);
          xlabel('Year');
          ylabel('Ring width (mm)');

          k=0;
          while k~=1;
            k=menu('Choose one: ','Zoom','Full','Laser','Menu');
            if k==1;
            elseif k==2
              figure('Color','k');
	      plot(yrx,x,'-',yry,y,'--');
	      title(['RW: Solid= ',Xfn,'   Dashed= ',Yfn]);
	      xlabel('Year');
	      ylabel('Ring Width (mm)');
            elseif k==3
	      eval('print -dps');
            elseif k==4
	      close
	      break
            end
          end ; % while loop on k
        end

      elseif k1==3;  %zoom analysis on transformed rw
        %  Now transform ring width and do same
        figure('Color','k');
        plot(yrw,w,'-',yrz,z,'--');
        title(['RW Change:  Solid= ',Xfn,'    Dashed = ',Yfn]);
        xlabel('Year');
        ylabel('Change in Ring Width (transformed)');
  
        k=1;
        while k~=4;
          [t,v]=ginput(2);  % get coordinates of end points of segment
          tbeg=ceil(t(1));
          tend=floor(t(2));
          if tbeg< max([yrw(1) yrz(1)])
            tbeg=max([yrw(1) yrz(1)]);
          end
          if tend > min([yrw(length(w))  yrz(length(z))]);
            tend=min([yrw(length(w))  yrz(length(z))]);
          end

          tseg=(tbeg:tend)';  % time segment for zoom
          iw = tseg-yrw(1)+1;  % index into w
          iz = tseg-yrz(1)+1;  % index into z

          % Kendall's tau and whether signif at 95%, 99%
          [T,tau]=kendtau(w(iw),z(iz));
          nT=length(w(iw));  % sample size for entering lookup table
          Tcrit=table1(cona12,nT); % get a row from the lookup table
          T95=Tcrit(2);  % 95 % signif level
          T99=Tcrit(4);  % 99% signif level
          if T>T99,
            sigT='**';
          elseif T>T95,
            sigT='*';
          else
            sigT=' ';
          end

          % Sign test
          [nagree,ntot,n95,n99]=signtest(w(iw),z(iz),1); 
          if n95=='Y' & n99=='N';
            sig='*';
          elseif n99=='Y';
            sig='**';
          else
            sig=' ';
          end

          % Binomial test
          [Nb,nsb,ps,pb]=binom1(w(iw),z(iz),1,.2);  % binomial test, 0.2 quantile
          if pb<=0.05 & pb >0.01;
            sigb='*';  % p-value 0.05 
          elseif pb <= 0.01;
            sigb='**';  % p-value 0.01
          else
            sigb=' ';
          end 
          figure('Color','k');
          plot(tseg,w(iw),'-',tseg,z(iz),'--');
          title(['RW Change:  Solid= ',Xfn,'    Dashed = ',Yfn]);
          xlabel ('Year');
          ylabel('Change in Ring width (transformed)');
          pltext(.5,.95,fonty,['tau = ',num2str(tau),' ',sigT]);
          pltext(.5,.90,fonty,['Sign: ',num2str(nagree),'/',num2str(ntot),' ',sig]);
          pltext(.5,.85,fonty,['Binom: ',num2str(nsb),'/',num2str(Nb),' ',sigb]);
          pltext(.1,.1,fonty,['N = ',int2str(length(w(iw)))]);

          k=0;
          while k~=1;
            k=menu('Choose one: ','Zoom','Full','Laser','Menu');
            if k==1;
            elseif k==2
              figure('Color','k');
	      plot(yrw,w,'-',yrz,z,'--');
	      title(['RW Change:  Solid= ',Xfn,'    Dashed = ',Yfn]);
	      xlabel('Year');
	      ylabel('RW change from previous year (scaled)')
            elseif k==3,
	      eval('print -dps');
            elseif k==4,
	      close
	      break
            end
          end;  % while loop on k
        end 

      elseif k1==4; % on outer while loop (k1)
      end
    end % of outer while loop  (k1)

    close all; % Close all graph windows.  

  end % End of outermost if loop

end  % End of outermost while loop


% End of file
