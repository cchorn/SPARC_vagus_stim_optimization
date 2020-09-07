function SI90 = classify_SI90_points(TS_data)

SI90 = {};
for i=1:length(TS_data)
    if TS_data(i)>=90
        SI90{i,1} = 'Above 90';
    else
        SI90{i,1} = 'Below 90';
    end
end

end