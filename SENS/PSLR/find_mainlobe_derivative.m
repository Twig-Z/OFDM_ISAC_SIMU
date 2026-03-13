%% 方案1：基于导数的主瓣边界检测
function mainlobe_indices = find_mainlobe_derivative(signal, peak_idx)
    N = length(signal);
    peak_val = signal(peak_idx);
    
    % 第一步：用较低阈值（如-6dB）粗略框定搜索范围
    coarse_threshold = peak_val * 0.5;  % -6dB
    
    left_search = peak_idx;
    while left_search > 1 && signal(left_search) > coarse_threshold
        left_search = left_search - 1;
    end
    
    right_search = peak_idx;
    while right_search < N && signal(right_search) > coarse_threshold
        right_search = right_search + 1;
    end
    
    % 第二步：在粗范围之外找第一个谷底作为真正边界
    % 向左找谷底（在left_search附近）
    left_idx = max(1, left_search - 1);
    while left_idx > 1
        if signal(left_idx) <= signal(left_idx-1) && ...
           signal(left_idx) <= signal(left_idx+1)
            break;
        end
        left_idx = left_idx - 1;
    end
    
    % 向右找谷底
    right_idx = min(N, right_search + 1);
    while right_idx < N
        if signal(right_idx) <= signal(right_idx-1) && ...
           signal(right_idx) <= signal(right_idx+1)
            break;
        end
        right_idx = right_idx + 1;
    end
    
    mainlobe_indices = left_idx:right_idx;
end