function [r] = rmse(data, estimate)
r=sqrt(sum((data(:)-estimate(:)).^2)/numel(data));
end

