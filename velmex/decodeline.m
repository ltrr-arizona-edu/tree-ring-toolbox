function meas = decodeline(tline)
%decodeline Turn a string from the stage controller into a measurement.
%   The string should contain the number at a fixed range of character
%   positions, followed by the units identifier "mm".

mmpos = findstr(tline, 'mm');
if isempty(mmpos)
    mmerr = MException('TreeRing:Measurement:NoMm', ...
        'Expecting the string \"mm\" in the line from the encoder.');
    throw(mmerr);
end;
meas = str2double(tline(mmpos-10:mmpos-1));
end
