% [S, cids] = parseKeggModel(fname, fmt)
%
% parses cellarray of reactions in KEGG format
%
% Arguments:
%   reactionStrings - a 1D cellarray containing strings of reactions
%   arrow - the string used as the 'arrow' in each reaction (default: '=')
% Return values:
%   S     - a stoichiometric matrix
%   cids  - the KEGG compound IDs in the same order as the rows of S
% Noor(2013)

function [S, rids, cids] = parseKeggModel(reactionStrings, arrow)

if nargin < 2
    arrow = '=';
end
if nargin < 3
    has_reaction_ids = false;
end

rids = {};
cids = {};
reactions = {};
for i = 1:length(reactionStrings)
    curr_line = reactionStrings{i};
    tokens = regexp(curr_line, '(\w+)\s+(.*)', 'tokens');
    rids = [rids,tokens{1}{1}];
    curr_line = tokens{1}{2};
    sprs = reaction2sparse(curr_line, arrow);
    cids = unique({cids{:}, sprs{:,1}});
    reactions = [reactions, {sprs}];
end

S = zeros(length(cids), length(reactions));
for i = 1:length(reactions)
    r = reactions{i};
    [~, idx] = ismember({r{:,1}}, cids);
    S(idx, i) = [r{:,2}];
end
