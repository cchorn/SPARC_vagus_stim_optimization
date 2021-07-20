function createnew_fig(cb, evendata)
% cb is the handle of the axes that was clicked. Click on the whitespace 
% within and axes and not on the line object and copy the axes object to 
% the new figure
hh = copyobj(cb,figure);
set(hh,'ButtonDownFcn',[]); %for the new figure assign the ButtonDownFcn to empty
set(hh, 'Position', get(0, 'DefaultAxesPosition')); %resize the axis to fill the figure
end