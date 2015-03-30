% jmd
% 8.15.12

%write cell array to file
%colnames is cell array
function fileID=writeCell(cell, filename, format, colnames, delimiter)

if nargin<3
    format=[repmat('%s ',1,size(cell,2)-1) '%s'];
end
if nargin<4
    colnames={};
end
if nargin<5
    delimiter=' ';
end

fileID = fopen(filename,'w');

if ~isempty(colnames)
    fprintf(fileID, [repmat('%s ',1,length(colnames)-1) '%s\n'], colnames{:}, delimiter);
end

for i=1:size(cell,1)
    fprintf(fileID, [format '\n'], cell{i,:}, delimiter);
end

fclose(fileID);