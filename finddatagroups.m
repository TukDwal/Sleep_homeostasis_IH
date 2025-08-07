function index = finddatagroups(data, value)

% Given a matrix numerical values, return the starting and ending index of
% every group of successive 'value'. If a group is constituted of one point, 
% starting and ending indexes are the same
% INPUT:
%   data: matrix of numerical values
%   value: value for which you are looking for groups
% OUTPUT: 
%   index: vector of length 2*n where n is the number of groups of values,
%   each pair of value is the starting and ending index of a group of
%   values
% ex:
% A = [0 0 1 1 1 0 0 1 1 0 1 0 0];
% out = finddatagroups(A, 1);
% out = [3 5 8 9 11 11]
%     ngroups = sum(diff(data==value));
    index = [];
    k1 = 1;
%     idx_group = 1;
    while k1 < length(data)
       if data(k1) == value
           k2 = k1+1;
           while data(k2) == value & k2 < length(data)
               k2 = k2+1;
           end
           index = [index k1 k2-1];
%            idx_group = idx_group + 1;
           k1 = k2;
       else
           k1 = k1 + 1;
       end
    end
end