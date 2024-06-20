
function idx_out = find_newline(str,idx_in)
% Find newline after certain index
idx_out = 1;
% while str(idx_in+idx_out) ~= newline
while str(idx_in+idx_out) ~= char(10)
    idx_out = idx_out+1;
end
end