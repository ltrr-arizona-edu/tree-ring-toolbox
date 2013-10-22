function [measuring, ix, yrthis] = pausemeasuring(yrorigin, yrthis, yrinmax, sport, reclen)
%pausemeasuring Prompt the user stop measuring, or change position.
%   Given the starting year, yrorigin, current year, yrthis, maximum
%   allowed year, yrinmax, serial port object, sport, and minimum input
%   record length, reclen, prompt the user to either stop measuring, or
%   switch to a different position in the series. The flag measuring
%   indicates whether or not measurement is to continue: if true, the user
%   has selected a new measurement position at year yrthis, with ix as the
%   updated index into the measurement data vector.

kmen4 = menu('Choose', 'Stop measuring', ...
    'Skip to a different part of series');
switch kmen4
    case 1
        disp('OK, Measurement stopped');
    case 2
        % skip to a differet part
        outofrange = true;
        while outofrange
            prompt={'Enter start year of segment to begin remeasuring):'};
            def = {int2str(yrthis + 2)};
            dlgTitle = 'Input year';
            lineNo = 1;
            answer = inputdlg(prompt, dlgTitle, lineNo, def);
            yrresume = str2double(answer{1});
            outofrange = (yrresume < yrorigin || yrresume > yrinmax);
            if outofrange
                uiwait(msgbox({'Year out of range; pick another!', ...
                    'If desired year out of range, increase maxlen or nyrlead', ...
                    'in calling function rwmeas.', ...
                    'Can change at opening prompt or the hardcoded defaults'}, 'Message', 'modal'));
            else
                yrthis = yrresume;
            end
        end
        reposition(sport, reclen);
end
ix = yrthis - yrorigin + 1;
measuring = (kmen4 == 2);
end

