function [xout, yrout, yrsagain] = measpart(xin, yrin, yrstart, yrsagain)

% yrsagain:  empty, or a rv of years flagged to require re-measuring; on return, empty if remeasured
%
%*** NOTES
%
% NaN rows in xin:  The input time series matrix generally has leading NaN rows to facilitate starting measurements before
% the start of existing measurements

close all;
clc;
% Set some year information
yrorigin = yrin(1);  % first year of the input time series matrix (tsm)
yrthis = yrstart; % start year for measurement
% Store the width components separately: total, early, late width
x1 = xin(:,1);
x2 = xin(:,2);
kwood = menu('Choose which part of first ring to begin measuring',...
    'Earlywood',...
    'Latewood');
sport = openserial();
% Initialize counters and while control
ix = yrthis - yrorigin + 1; % index into x for storing data for yrthis
reclen = 15; % input string length
measuring = true;
dnum = 0;
reposition(sport, reclen);
% Measurement while loop
while measuring
    try
        count = 0;
        while count < reclen
            [g, count] = fgets(sport); % store measurement (a string) in g
        end
        prevnum = dnum;
        dnum = decodeline(g);
        % measurement is difference of cumulative reading
        w = dnum - prevnum;
        if w < 0
            % negative delta -- meaning you reversed crank on velmex
            dnum = 0;
            ix = ix - 1;
            disp(num2str([yrthis x1(ix) x2(ix) x1(ix)+x2(ix)]));
            [measuring, ix, yrthis] = pausemeasuring(yrorigin, yrthis, max(yrin), sport, reclen);
            if measuring
                kwood = menu('Resume measuring at ', 'Earlywood', 'Latewood');
                pause(1)
            end
        else
            % measurement is non-negative
            [x1, x2, ix, kwood, yrthis, yrsagain] = measuretwo(x1, x2, ix, w, kwood, yrorigin, yrthis, yrsagain);
        end
    catch EX
        disp(EX.message);
        closeserial(sport);
        rethrow(EX);
    end
end;
[xout, yrout] = tidyouttwo(x1, x2, yrin);
closeserial(sport);
end

function [x1, x2, ix, kwood, yrthis, yrsagain] = measuretwo(x1, x2, ix, w, kwood, yrorigin, yrthis, yrsagain)
% Measure earlywood and latewood widths.
if kwood == 1;
    x1(ix) = w;
    kwood = 2;
else
    x2(ix) = w;
    kwood = 1;
    disp(num2str([yrthis x1(ix) x2(ix) x1(ix)+x2(ix)]));
    if rem(yrthis, 10) == 0
        if ispc  % Test OS type to play sound on decade.
            beep;
        elseif isunix
            !aplay -q /usr/local/MATLAB/MLB/velmex/decade-beep.wav &
        end
    end
    % Clear remeasure flag if year was flagged for re-measurement
    if ~isempty(yrsagain)
        yrsagain(yrthis == yrsagain) = [];
        if isempty(yrsagain)
            yrsagain = [];
        end
    end
    yrthis = yrthis + 1;
    ix = yrthis - yrorigin + 1;
end
end

function [xout, yrout] = tidyouttwo(x1, x2, yrin)
% Strip leading and trailing NaN values and reorganize output arrays.
x = x1 + x2;
X = [x1 x2 x];
mx = size(X,1);
% Strip trailing NaN
L = (all(isnan(X')))';
h = ones(mx, 1);
h(L) = NaN;
h = trailnan(h);
mh = length(h);
X = X(1:mh,:);
yrin = yrin(1:mh);
% strip leading NaN
h = flipud(h);
h = trailnan(h);
mh = length(h);
X = flipud(X);
X = X(1:mh,:);
X = flipud(X);
yrin = flipud(yrin);
yrin = yrin(1:mh);
yrin = flipud(yrin);
xout = X;
yrout = yrin;
end
