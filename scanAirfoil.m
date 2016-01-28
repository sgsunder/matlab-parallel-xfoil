function airfoil = scanAirfoil(filename)
% Open File
fID = fopen(filename);
% Read Header, Get Name
header      = fgetl(fID);
has_header  = cellfun(@isempty, ...
    textscan(header, '%f%f', 'MultipleDelimsAsOne', true));
has_header  = has_header(1);
if has_header
    airfoil.name = strtrim(header);
else
    [~, airfoil.name, ~] = fileparts(filename);
end
% Read Coordinate Data
frewind(fID);
if ~isnumeric(has_header)
    has_header = 1;
end
data = textscan(fID, '%f%f', 'HeaderLines', has_header, ...
    'MultipleDelimsAsOne', true);
if isempty(data{1}) || isempty(data{2})
    data = textscan(fID, '%f,%f', 'HeaderLines', has_header, ...
        'MultipleDelimsAsOne', true);
end
fclose(fID);
% Separate Upper and Lower Coordinates
X = data{1};
Y = data{2};
% Normalize lengths and center leading edge
scale = 1/(max(X) - min(X));
X = (X - min(X)) .* scale;
Y = Y .* scale;
% Split into Upper and Lower Sections
iLE = find(X==0, 1, 'first');
airfoil.UX = X(iLE:-1:1);
airfoil.UY = Y(iLE:-1:1);
airfoil.LX = X(iLE:end);
airfoil.LY = Y(iLE:end);
end

