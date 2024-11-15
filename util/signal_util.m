1; % This is not a function file

function [ts, as] = generate_signal(tt, yy, dt)
    % ts: time points
    % aa: altitude
    % tt: time key points
    % yy: altitude key points
    ts = min(tt):dt:max(tt);
    as = interp1(tt, yy, ts, 'pchip');
end

function [yy] = merge_signal(yy, tt2, yy2, tt2_offset, dt)
    % signal 1 - base signal, signal 2 - signal to merge.
    % Will approximate signal 2, and merge
    % tt2 - signal 2's base points
    tt2_pos = (tt2_offset) / dt;
    [tt2, yy2] = generate_signal(tt2, yy2, dt);
    yy(tt2_pos:tt2_pos + numel(yy2) - 1) += yy2;
end

function val = clamp(low, value, high)
    val = min(high, max(low, value));
end
