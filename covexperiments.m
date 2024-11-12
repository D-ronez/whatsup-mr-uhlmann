function [] = main()
    source("signal_util.m")
    % Trajectory base points
    tt = [0.0, 14.280601288475305, 22.440944881889763, 33.321403006442374, 35.36148890479599, 102.68432355046528, 106.76449534717251, 138.04581245526128, 156.4065855404438, 173.40730136005726, 184.28775948460986, 194.48818897637796, 207.40873299928418, 213.52899069434503, 218.96921975662133, 224.40944881889763, 231.88976377952756, 245.49033643521832, 254.33070866141733, 265.2111667859699, 275.411596277738, 285.6120257695061, 292.41231209735145, 295.8124552612742, 310.09305654974946, 318.93342877594847, 320.97351467430207, 361.0952040085898, 371.2956335003579, 388.29634931997134, 427.05798138869005, 445.4187544738726, 452.21904080171794, 465.8196134574087, 473.9799570508232, 486.2204724409449, 492.34073013600573, 501.86113099498925, 517.5017895490337, 532.4624194702935, 541.982820329277, 552.8632784538296, 564.4237652111668, 573.2641374373658, 584.8246241947029, 596.3851109520401, 637.8668575518969, 659.6277738010021, 722.8704366499642, 728.3106657122405, 765.0322118826056, 781.3528990694344, 792.9133858267717, 797.6735862562634, 846.6356478167502, 852.0758768790265, 860.2362204724409, 870.436649964209, 877.2369362920543, 890.8375089477452, 928.9191123836792];
    yy = [41.58536585365854, 47.92682926829268, 52.073170731707314, 57.19512195121951, 58.90243902439024, 58.90243902439024, 52.4390243902439, 50.48780487804878, 51.34146341463415, 51.829268292682926, 51.70731707317073, 53.048780487804876, 53.90243902439024, 55.48780487804878, 57.80487804878049, 57.92682926829268, 56.34146341463415, 57.80487804878049, 56.34146341463415, 57.92682926829268, 53.90243902439024, 52.92682926829268, 53.292682926829265, 51.951219512195124, 51.70731707317073, 51.46341463414634, 50.609756097560975, 50.36585365853659, 51.46341463414634, 50.853658536585364, 49.8780487804878, 51.34146341463415, 52.3170731707317, 52.073170731707314, 53.292682926829265, 53.90243902439024, 57.073170731707314, 56.34146341463415, 58.65853658536585, 56.46341463414634, 56.951219512195124, 56.829268292682926, 53.53658536585366, 52.5609756097561, 52.3170731707317, 51.21951219512195, 50.48780487804878, 51.70731707317073, 51.58536585365854, 50.0, 49.8780487804878, 51.09756097560975, 53.170731707317074, 60.0, 59.63414634146341, 54.63414634146341, 54.63414634146341, 44.390243902439025, 44.26829268292683, 41.58536585365854, 41.58536585365854];
    [tt, yy] = generate_signal(tt, yy, 0.5);
    dt = tt(2) - tt(1)

    % GNSS
    % Add noise
    gsigma = 0.3;
    gyy = yy + randn(size(yy)) * 0.1;
    % Add sag
    tsag = [0 5 10 20];
    ysag = [0 -8 5 0];
    gyy = merge_signal(gyy, tsag, ysag, 2 * 60, dt);
    % Add another sag
    tsag = [0 2 4 8];
    ysag = [0 -8 5 0];
    gyy = merge_signal(gyy, tsag, ysag, 30, dt);
    % Add another sag
    tsag = [0 2 3 5];
    ysag = [0 -1 1 0];
    gyy = merge_signal(gyy, tsag, ysag, 600, dt);

    % Baro
    % Add noise
    byy = yy + randn(size(yy)) * 1.3;
    % Add drift, 5m
    byy = byy + linspace(0, 5, numel(yy));

    % Stats
    [cyy, corr, vyyg, vyyb] = get_cov_window(tt, gyy, byy, 20);

    % Filtering
    sigma_frac = vyyg ./ vyyb;
    % Trust boost
    boost = -sigma_frac + 0.2;
    gnss_score = (boost + vyyb) ./ (vyyb + vyyg);
    for i = 1:numel(gnss_score)
        gnss_score(i) = clamp(0, gnss_score(i), 1);
    end
    baro_score = ones(size(vyyb)) - gnss_score;
    fyy = gnss_score .* gyy + baro_score .* byy;

    ax1 = subplot(2, 1, 1)
    hold on
    grid on
    plot(tt, gyy, '-')
    plot(tt, byy, '-')
    plot(tt, yy, '-', 'LineWidth', 2)
    plot(tt, fyy, '-', 'LineWidth', 2)
    legend("GNSS", "Baro", "Ground truth", "Filtered")
    ax2 = subplot(2, 1, 2)
    hold on
    grid on
    plot(tt, cyy)
    plot(tt, sigma_frac)
    legend("Cov. b/w measures", "sigma GNSS / sigma BARO")
    linkaxes([ax1, ax2], "x")
    % subplot(3, 1, 3)
    % plot(tt, corr)
    % hold off
    % legend("Correlation")
end

function [cyy, corr, vyy1, vyy2] = get_cov_window(tt, yy1, yy2, window_size)
    % Builds sliding window covariance
    % yy1 - signal 1
    % yy2 - signal 2
    % tt - time
    % cyy - sliding window covariance
    % corr - correlation (window)
    % vyy1 - variance of signal 1

    window = ones(size(yy1)) * window_size;
    window(1:window_size) = linspace(1, window_size, window_size);
    cyy = zeros(size(yy1));
    vyy1 = zeros(size(yy1));
    vyy12 = zeros(size(yy1));
    corr = zeros(size(yy1));
    for i = 2:numel(yy1)
        win = window(i);
        % Covariance
        cyy(i) = cov(yy1(i - win + 1:i), yy2(i - win + 1:i));
        % Standard deviations
        sigma1 = std(yy1(i - win + 1:i));
        sigma2 = std(yy2(i - win + 1:i));
        vyy1(i) = sigma1;
        vyy2(i) = sigma2;
        % Correlation
        corr(i) = cyy(i) / (sigma1 * sigma2);
    end
end
