function obj = openserial()
%openserial Open the serial port.
%   The settings defined here must match those on the stage controller.

% Just in case, make sure serial port not assigned to some object

scheck = instrfind;
if ~isempty(scheck)
    fclose(scheck);
    delete(scheck);
    clear scheck;
end
% Create serial port object  
obj = serial('/dev/ttyS0'); % For Ubuntu systems.
%obj = serial('com1');  % Uncomment this and comment the above line for win32.
% Set parameters and open serial port
set(obj, ...
    'Parity', 'none', ...
    'StopBits', 1, ...
    'DataBits', 8, ...
    'BaudRate', 9600', ...
    'terminator', 'CR/LF', ...
    'TimeOut', 600);
fopen(obj);

end

