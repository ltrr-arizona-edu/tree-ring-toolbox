function closeserial(obj)
%closeserial Close a serial port.
%   Given obj, a serial port object, this closes the port and destroys the
%   object.

fclose(obj);
delete(obj);
clear obj;
end

