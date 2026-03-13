%% 子函数: 可视化
function visualize_PSLR(tau_vec, fd_vec, chi, chi_0delay, chi_0doppler, ...
                        mainlobe_idx_d, mainlobe_idx_f, ...
                        peak_idx_d, peak_idx_f, ...
                        sidelobe_idx_d, sidelobe_idx_f)
figure('Position', [100 100 1400 500]);
subplot(1,2,1);
mesh(tau_vec, fd_vec, chi');
xlabel('时延 τ (s)'); ylabel('多普勒 f_d (Hz)'); zlabel('幅度');
title('模糊函数χ(τ, f_d) 立体');
axis tight; view(45, 30);
    
    
    
    % 完整模糊函数
    subplot(1,2,2);
    imagesc(fd_vec, tau_vec, chi );
    set(gca, 'YDir', 'normal');
    xlabel('多普勒 f_d (Hz)');
    ylabel('时延 τ (s)');
    title('模糊函数 χ(τ, f_d) 平面');
    colorbar;
    %clim([-60 0]);
    

    figure('Position', [100 100 1400 500]);
    % Zero-Doppler切片
    subplot(1,2,1);
    plot(tau_vec, chi_0doppler + 1e-12, 'b-', 'LineWidth', 1.5);
    hold on;
    % 标记主瓣
    plot(tau_vec(mainlobe_idx_d), ...
         chi_0doppler(mainlobe_idx_d) + 1e-12, 'g.', 'MarkerSize', 15);
    % 标记峰值
    plot(tau_vec(peak_idx_d), ...
         chi_0doppler(peak_idx_d) + 1e-12, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    % 标记最高旁瓣
    plot(tau_vec(sidelobe_idx_d), ...
         chi_0doppler(sidelobe_idx_d) + 1e-12, 'mx', 'MarkerSize', 12, 'LineWidth', 2);
    
    xlabel('时延 τ (s)');
    ylabel('幅度');
    title('Zero-Doppler切片 χ(τ, 0)');
    legend('切片', '主瓣', '峰值', '最高旁瓣', 'Location', 'best');
    grid on;
    %ylim([-60 5]);
    
    % Zero-Delay切片
    subplot(1,2,2);
    plot(fd_vec, chi_0delay + 1e-12, 'b-', 'LineWidth', 1.5);
    hold on;
    % 标记主瓣
    plot(fd_vec(mainlobe_idx_f), ...
         chi_0delay(mainlobe_idx_f) + 1e-12, 'g.', 'MarkerSize', 15);
    % 标记峰值
    plot(fd_vec(peak_idx_f), ...
         chi_0delay(peak_idx_f) + 1e-12, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    % 标记最高旁瓣
    plot(fd_vec(sidelobe_idx_f), ...
         chi_0delay(sidelobe_idx_f) + 1e-12, 'mx', 'MarkerSize', 12, 'LineWidth', 2);
    
    xlabel('多普勒 f_d (Hz)');
    ylabel('幅度');
    title('Zero-Delay切片 χ(0, f_d)');
    legend('切片', '主瓣', '峰值', '最高旁瓣', 'Location', 'best');
    grid on;
    %ylim([-60 5]);
end