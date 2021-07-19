function [row,column] = gridLoc(grid,chanNum)
%gridLoc finds outputs the row and column of a channel when given a channel
%number and a grid outline
%   
    temp_loc = chanNum == grid;
    [row, column] = find(temp_loc);
end

