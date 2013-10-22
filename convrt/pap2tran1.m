function pap2tran1
% pap2tran:  reformat figure from paper to transparency

prompt={'Enter scale factor for font sizes:',...
        'Enter scale factor for line widths:'};
   def={'1.2','3'};
   dlgTitle='Expansion Factors';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
   factfs1 = str2num(answer{1});
   factlw1 = str2num(answer{2});

   

   
h=findobj; % get all object handles

nh=size(h,1);

for n = 1:nh;
    hthis=h(n);
    htype = get(hthis,'Type');
    htag = get(hthis,'Tag');
    
    if strcmp(htype,'text');
        fs = get(hthis,'FontSize');
        set(hthis,'FontSize',fs*factfs1);
        set(hthis,'FontName','arial');
    end;
    
    if strcmp(htype,'line');
        lw = get(hthis,'LineWidth');
        set(hthis,'LineWidth',lw*factlw1);
    end;
    
    
end;

% Find out how many axes
haxe=get(gcf,'Children');
nkids=length(haxe);
% Loop over children, and beat axes only
for n = 1:nkids
    hthis=haxe(n);
    if strcmp('axes',get(hthis,'Type'))
        % Tick labels
        fsaxes=get(hthis,'FontSize');
        set(hthis,'FontSize',factfs1*fsaxes);
        set(hthis,'FontName','Arial');
        
        % labels of axes
        hlabely=get(hthis,'YLabel');
        hlabelx=get(hthis,'XLabel');
        fslabel=get(get(hthis,'YLabel'),'FontSize');
        set(hlabely,'FontSize',factfs1*fslabel);
        set(hlabelx,'FontSize',factfs1*fslabel);
        set(hlabely,'FontName','Arial');
        set(hlabelx,'FontName','Arial');
        
        % title
        htitle=get(hthis,'Title');
        fstitle=get(htitle,'FontSize');
        set(htitle,'Fontsize',fstitle*factfs1);
    end;
    
end; % 