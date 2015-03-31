% jmd
% 6.23.12
% sbmlCode2ascii.m

% need to regexp left-to-right, to match '__' w/ correct place in name
% maybe can use compList in readCbModel to avoid '__93_[_]' ending
% need to also deal w/ _\d at beginning of some names

% sbmlCell is a cell array of strings
function ascii=compliantSbmlNames2ascii(sbmlCell)
% translates a string from a compliant SBML that has gone thru readCbModel

% in sbml, names can't start w/ a number, so there's an extra '_' in front
sbmlCell=regexprep(sbmlCell, '^_(\d)', '$1');
% readCbModel assumes end has [.], where '.' is a one-letter abbreviation
% of the compartment name, so we need to undo the effects of this
sbmlCell=regexprep(sbmlCell, '\[(.)\]$', '_$1');

ascii=cell(length(sbmlCell),1);
for i=1:length(sbmlCell)
    % non-SBML characters have been replaced by '__decimal ascii code__'
    % this code uses 2-to-3 digits to represent printable non-SBML characters
    [matchstart,matchend,tokenindices,matchstring,tokenstring,tokenname, ...
            splitstring] = regexpi(sbmlCell{i}, '__(\d{2,3})__');
    % splitstring holds parts of input string delimited by expression,
    % so res takes its substrings and leaves room for converted codes
    res=cell(1, 2*length(splitstring)-1);
    res(1:2:length(res))=splitstring;
    % loop through & transform ascii code to ascii using str2double
    for t=1:length(tokenstring)
        res{2*t}=char(str2double(char(tokenstring{t})));
    end
    ascii{i,1}=[res{:}];
end