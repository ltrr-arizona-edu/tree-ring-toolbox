function coreplt1(X,nms,yrs,S);
% coreplt1;Subfunction to coreview.m. Plots ring widths, age trend, index, for 1 core
% CALL: coreplt1(X,nms,yrs,S);
% 
% 
% Choose core
kcore = menu('Choose core',cellstr(nms));

% Pull ringwidth series
yrgo1 = yrs(kcore,1);  yrsp1 = yrs(kcore,2);  % start, end year for ring width
nyrs = yrsp1-yrgo1+1; % number of years of ringwidt
nmcore =char(nms(kcore)); % core name
igo1= yrs(kcore,3); % start row of ring-width in X
x1 = X(igo1:igo1+nyrs-1); % ring width
yr1 = (yrs(kcore,1):yrs(kcore,2))'; % year vector for ring-width

% Plot ring-width
figure(1);
hp1 = plot(yr1,x1);
return;


