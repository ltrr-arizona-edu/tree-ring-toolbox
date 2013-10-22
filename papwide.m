% papwide.m  size a full-page figure so that width exactly fills paper width
%
% You control by width because might need greater length of paper than heighth
% of figure to accomodate a figure caption, say, in word

%-- Prompt for optional margin at left and right of figure
prompt={'Enter margin (in):'};
def={'0.5'};
dlgTitle='Margin to keep at L and R of figure';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
sidemarg=str2num(answer{1});

%-- Prompt for optional margin at top and bottom of figure
prompt={'Enter margin (in):'};
def={'0.5'};
dlgTitle='Margin to keep at T and B of figure';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
tbmarg=str2num(answer{1});


%-- Prompt for possible room at bottom of figure for caption
kmen3 = menu('Choose one',...
   'Figure will have a caption',...
   'Figure will not have a caption');
if kmen3==1;
   prompt={'Enter heigth (in):'};
   def={'1'};
   dlgTitle='Height of Area for Optional Caption';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
   dcaption=str2num(answer{1});
else;
   dcaption=0;
end;

%--- Menu to find out paper type
kmen1 = menu('Choose one','Letter','Legal','Super A3/B');
if kmen1==1;
   psize=[8.50 11]; % paper size letter
elseif kmen1==2;
   psize=[8.5 14]; % legal
elseif kmen1==3;
   psize=[12.96 19.02];
end;

%-- Figure whether to use portrait or landscape
figpos=get(gcf,'Position'); % Position of figure as created
WLratiofig = figpos(3)/figpos(4);

%-- Menu to allow you to specify landscape or portrait
kmen2 = menu('Choose one',...
   'Landscape',...
   'Portrait',...
   'Let width/height ratio determine orientation');
if kmen2==1; 
   set(gcf,'PaperOrientation','Landscape');
elseif kmen2==2;
   set(gcf,'PaperOrientation','Portrait');
else; % If figure wider than high let width fill landscape, otherwidse portrait
   if WLratiofig>1;
      set(gcf,'PaperOrientation','Landscape');
   else;
      set(gcf,'PaperOrientation','Portrait');
   end;
end;

if strcmp('landscape',get(gcf,'PaperOrientation'));
   psize=fliplr(psize); % psize is now [width height] of paper, adjusted for orientation
end;

% Compute paper position.  This will be limited either by the figure or the paper
rat1=(psize(1)-2*sidemarg)/(psize(2)-dcaption-2*tbmarg)  ;
if rat1 > WLratiofig; % let height of paper limit figure size
   p4 = psize(2)-dcaption-2*tbmarg; % heigth of printed figure
   p3 = p4 * WLratiofig; % width of printed figure
else;
   p3 = psize(1) - 2*sidemarg;
   p4 = p3 / WLratiofig;
end;

p1 = sidemarg;
p2 = tbmarg;

set(gcf,'PaperPosition',[p1 p2 p3 p4]);
