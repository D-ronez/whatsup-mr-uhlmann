function [] = main()
    % Trajectory base points
    tt = [0 5 11 11.5 12 12.5 13 14] * 60;
    yy = [0 50 40 50 55 51 10 0];
    [tt, yy] = generate_signal(tt, yy, 0.5);
    dt = tt(2) - tt(1)

    % GNSS
    % Add noise
    gsigma = 0.3;
    gyy = yy + randn(size(yy)) * 0.1;
    % Add sag
    tsag = [0 5 10 15]
    ysag = [0 -10 5 0]
    gyy = merge_signal(gyy, tsag, ysag, 2 * 60, dt)

    % Baro
    % Add noise
    byy = yy + randn(size(yy)) * 0.4;
    % Add drift, 5m
    byy = byy + linspace(0, 5, numel(yy));

    hold on
    plot(tt, gyy, '--')
    plot(tt, byy, '--')
    plot(tt, yy, '-', 'LineWidth', 2)
    legend("GNSS", "Baro", "Ground truth")
    hold off
end

function [tt, yy] = add_noise(tt, yy, var)
end

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
