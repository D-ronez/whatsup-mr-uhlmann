function [series] = cooldown_cumsum(series, cooldown)
	series(1) = cooldown * series(1);
	for i = 2:numel(series)
		series(i) = cooldown * (series(i) + series(i - 1));
	end
end
