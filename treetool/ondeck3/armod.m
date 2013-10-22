function [EX,P] = armod(IX,S,yrs,nms,nnmx)
% armod: fit AR models to core ring width, under control of corei.m
% [EX,P] = armod(IX,S,yrs,nms,nnmx);
% Last revised 9-2-99
%
% Utility function to fit AR models to core indexh series and prewhiten them. Dedicated 
% function used by corei.m.  Returns AR model parameters and stats as well as residuals. 
% Appends information to a .mat tree-ring storage file.
%
%*** INPUT
% 
% IX (? x 1)r sov of standard core indices
% S (? x ?)r conventional core info in tree-ring .mat files
% yrs (? x 3)i  years and row-index data
% nms (? x ?)s core ids
% nnmx (1 x 1)i maximum allowable AR order
%
%*** OUTPUT
%
% EX (? x 1)r residual core indices (same	length as IX)
% P (? x 16) - Matrix containing AR parameters and other statistics
%
%*** REFERENCES -- NONE
%*** UW FUNCTIONS CALLED
% resid1
% acf
% pacf
%
%*** TOOLBOXES NEEDED
% system identification

[ns,dum]=size(nms);

% Find the total # of blocked segments
[rs cs]=size(S);
dsm=0;
for j=1:rs,
  dsm=dsm+length(find(S(j,cs-9:cs)));
end
tblkn=dsm/2;
   
% Initialize the P matrix
rpn=ns+tblkn;
P=ones(rpn,16)*NaN;
pidx=0;

% Initialize EX
EX=ones(length(IX),1)*NaN;

txtmsg2=['You have specified max AR order to consider as ' int2str(nnmx)];
hmsg2 = msgbox(txtmsg2);
pause(1);
close(hmsg2);

% Prompt if user wants to view the graph
klook=questdlg('WANT TO VIEW PLOTS OF ORIGINAL SERIES AND AR RESIDUALS AS MODELING PROCEEDS?');

switch klook;
case 'No';
   % no action
case 'Cancel';
   klook='No';
case 'Yes';
   hmsg=msgbox('Press Mouse Button to Move from Plot to Next Plot'); 
   pause(3);
   close(hmsg);
end

if ~isequal(klook,'Yes');
   hmsg=msgbox('Hang on... doing the AR modeling');
   pause(1);
   close(hmsg);
end

% Loop through for P and ex
for scid=1:ns,
   if S(scid,1)~=0,
      yrv=(yrs(scid,1):yrs(scid,2))';
      xv=IX(yrs(scid,3):yrs(scid,3)+yrs(scid,2)-yrs(scid,1));
      lzn=isnan(xv);
      zn=xv;
      zn(lzn)=[];
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
      
      for k=1:2:length(lz)-1,
         xvt=xv(lz(k)+1:lz(k+1));
         pidx=pidx+1;
         if length(xvt)>=20,
            P(pidx,1)=scid;
            P(pidx,2)=yrv(lz(k)+1);
            P(pidx,3)=yrv(lz(k+1));
            
            % Selecting the model order
            xvtt=detrend(xvt,0);	      % Remove the mean from the series
            v=arxstruc(xvtt,xvtt,[1:nnmx]');
            nn=selstruc(v,'aic');
            
            P(pidx,4)=nn;
            
            th=ar(xvtt,nn);
            temex=resid1(xvtt,th);
            
            %P(pidx,5)=var(temex)/var(xvtt);
            varrat=th(1,1)/var(xvtt); % ratio of variance of residuals to original variance
            P(pidx,5)=varrat;
            P(pidx,6:5+nn)=th(3,1:nn);
            
            temex=temex+mean(zn)-mean(temex);
            temex(1:nn)=ones(nn,1)*NaN;
            EX(yrs(scid,3)+lz(k):yrs(scid,3)+lz(k+1)-1)=temex;
            
            [r,SE2,r95]=acf(zn,3);
            P(pidx,11:13)=r;
            [phi,SE2] = pacf(zn,3);
            P(pidx,14:16)=phi;
            
            if isequal(klook,'Yes'); % want to view plots
               % Plot the rw data with scid core ID
               %hf0=figure('Color','w','Units','normal');
               hf0=figure('Units','normal');
               yrvt=yrv(lz(k)+1:lz(k+1));
               hp1=plot(yrvt,xvt,'b',...
                  yrvt,temex,'m--');title(['AR Modeling of ',nms(scid,:)]);
               %set(hp1(2),'Linewidth',1.5);
               ht1=text('Position',[0.1 0.9],'String',['AR order = ',num2str(nn)],...
                  'Unit','Normalized');
               ht4=text('Position',[0.1 0.85],'String',...
                  ['Var(resids)/Var(original) = ' sprintf('%4.2f',varrat)],...
                  'Unit','Normalized');
               ht2=text('Position',[0.6 0.9],'String','Dashed (magenta) - Residual','Unit','Normalized');
               ht3=text('Position',[0.6 0.85],'String','Solid (blue) - Original','Unit','Normalized');
               xlabel('Year');ylabel('Index');
               %jd=jdisp('Hit return to continue',[0.05 0.3 0.7 0.1]);
               %pause;
               %close(jd);
               T=waitforbuttonpress;
               close(hf0);
            end
         else
            error(['This segment of ',nms(scid,:),' is shorter than 20. Can''t ',...
                  'process. NaN''s will appear in the P matrix']);
         end
      end
   else
      pidx=pidx+1;
   end		% End of if loop
   disp(['   ' nms(scid,:) ' modeled as AR(' int2str(P(scid,4)) ')']);
end
close all



% End of file
