function SNR_dB=constellation_evm(data_aftereq,limit)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
data_aftereq_real=real(data_aftereq);
data_aftereq_imag=imag(data_aftereq);
figure(1)
subplot(2,3,5)
hold off
scatter(data_aftereq_real,data_aftereq_imag,5,'filled');
hold on 
modulation_x=[-sqrt(2)/2,-sqrt(2)/2,sqrt(2)/2,sqrt(2)/2];
modulation_y=[-sqrt(2)/2,sqrt(2)/2,-sqrt(2)/2,sqrt(2)/2];
scatter(modulation_x,modulation_y,20,'filled','red');

% parfor idx=1:length(data_aftereq)
%     for constellation_point=1:4
%         delta_x=(idx)=data_aftereq_real(idx)-modulation_x(constellation_point);
%         delta_y=(idx)=data_aftereq_imag(idx)-modulation_y(constellation_point);
%         delta_x^2+
%     end
% end

difference_all=zeros(4,length(data_aftereq));
difference_all(1,:)=abs(data_aftereq-(sqrt(2)/2+1j*sqrt(2)/2)*ones(1,length(data_aftereq)));
difference_all(2,:)=abs(data_aftereq-(sqrt(2)/2-1j*sqrt(2)/2)*ones(1,length(data_aftereq)));
difference_all(3,:)=abs(data_aftereq-(-sqrt(2)/2+1j*sqrt(2)/2)*ones(1,length(data_aftereq)));
difference_all(4,:)=abs(data_aftereq-(-sqrt(2)/2-1j*sqrt(2)/2)*ones(1,length(data_aftereq)));
difference_min=min(difference_all);
Power_error=mean(difference_min.^2);
Power_ref=1;
EVM_linear=sqrt(Power_error/Power_ref);
EVM_dB=20*log10(EVM_linear);
% SNR_linear=(1)/sqrt(EVM_linear);
PAPR=0;
SNR_dB=-EVM_dB+3+0;

xlim ([-limit limit]);
ylim ([-limit limit]);
xlabel('I 路','FontWeight','bold','FontName','fangsong','FontSize',15);
ylabel('Q 路','FontWeight','bold','FontName','fangsong','FontSize',15);
title('星座图','FontWeight','bold','FontName','fangsong','FontSize',15);
end

