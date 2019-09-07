function [r] = mse(data, estimate)
error2 = (data - estimate).^2;
r=mean(error2,3);
end

