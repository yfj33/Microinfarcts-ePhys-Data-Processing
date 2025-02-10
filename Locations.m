clc;clear;close all
Distance=[423 493 1147 1813 640 380 1183 507 787 1186 490 0 303 906 1183 506 786 1186];
f=figure
histogram(Distance,20);
ax=gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
%saveas(f,)
%less than 600 um as close, list recording data
R_far=[nan 600.6618  563.0747  539.8323  508.5250  498.5672  499.2388  525.6171  558.5511 nan 485.9996  498.9063  487.1097  503.6580  486.4321  549.9910  436.3614  460.8612  368.0310;
;nan 20.9278   38.0107   34.0639   30.7575   26.9073   21.9065   22.2680   18.0290 22 23.5560   40.2971   46.2743   34.5186   54.9666   78.0704   30.3715   20.3450   54.3186;
nan   56.7446   61.0296   52.1261   74.8647   83.5205   62.8421   97.7951   85.3201 nan  53.0778   94.6004   87.3896   52.6003  109.4976   83.6401   55.5754   25.8509   66.5910;
 216.7654  137.4781  145.2320  109.5984  308.9639  146.9179  158.3823  135.6809  173.6347 190 214.5000  166.6158  146.7815  143.8789  244.2297  205.0709  244.6926  380.0924  242.9528];

%600-1000um as middle 
R_middle=[nan 194.0892  236.0785  306.4496  286.3668  301.5913  284.5451  247.1267  216.9664  223.1903  189.1721  116.4634  160.5643  213.6922  179.4243  226.4642  214.0657 nan nan; 
     250.1490  208.5556  208.0417  216.8856  216.4248  249.1278  210.6551  188.0126  208.2898  nan 243.2586  204.1016  197.5067  228.3312  217.7000  178.3652  214.2836  198.6962  198.2976;
      169.6422  173.9789  281.9024  149.9920  172.3770  173.4538  151.9290  167.4876  149.2422  130.7835  177.7677  154.0343  202.3710  164.7373 nan nan nan nan nan;
     234.9220  257.3258  225.9462  235.4195  240.3313  254.9110  224.1377  253.7788  236.9756  258.5732  301.7829  261.5201  305.2438  215.7574 nan nan nan nan nan; 
      167.1261  188.0306  344.8623  189.7180  199.3730  208.2072  134.5434  162.9442  140.1275  146.6293  181.9222  163.3212  192.4753  162.3844 nan nan nan nan nan;
];

%larger than 1000 um as far 
R_close=[nan 80.8718   46.7920   65.4386   53.7072   34.5517   23.0906   44.2359   31.2509  nan 34.7094   65.4987   62.5807   43.0183   85.4403   35.4517   55.4671   83.0756  101.5410
;nan 422.3456  351.2126  259.5922  358.3990  232.2671  178.7126  156.5832  207.7051 nan 247.5733  248.6869  375.8765  277.8034  391.8476  367.0662  349.5361  284.1922  321.6347;
 nan 340.3769  399.7500  428.2563  351.9121  164.4409  219.5171  211.3725  248.9760  379.7858  384.4937  357.5281  405.4478  437.7055  382.2871  341.3896  347.9611 nan nan;
  nan 377.7278 428.4740 475.6334 4.4690 1.6661 0 3.5557 0 1.0383 4.4714 51.6639 59.0303 76.5425 49.5418 77.3204 84.8511 nan nan;    
  nan 416.3790  499.1955  545.7645  142.6540   67.0661   42.8727   46.8288   78.7222  166.6529  166.7465  187.2095  193.2308  233.7511  164.5414  229.9179  208.5189 nan nan;
  172.5262  153.8206  166.6381  161.3128  143.8792  134.7816  161.8131  153.7689 153.9136 nan 142.2353  176.1625  165.8979  168.2860  208.6638  177.3609  162.7912  165.0206  184.2287;
325.0828  377.8928  313.5868  324.8822  329.4891  366.3898  247.3900  290.9391  304.7969  310.4047  399.2587  439.6382  418.1765  372.7419 nan nan nan nan nan;
 108.5794  108.4078   77.0918   81.6454    7.1153    4.2143    1.4476    3.3403    3.2400   11.4536    5.3941    8.4947    5.5263    3.7000    5.2396    6.1330    2.5664    7.9422    9.7100;
  290.4111  276.2624  232.9827  248.7020  129.6022  117.2338  100.1626  132.4250  156.4976  183.4150  176.9015  202.1676  258.3334  250.6222  234.8284  303.7333  300.8256  348.1053  313.6260];


for i=1:size(R_far)
    R_far_norm(i,:)= R_far(i,:)/nanmean(R_far(i,1:4));
   
end

for i=1:size(R_middle)
   % R_far_norm(i,:)= R_far(i,:)/nanmean(R_far(i,1:4));
     R_middle_norm(i,:)=R_middle(i,:)/nanmean(R_middle(i,1:4));
    %R_close_norm(i,:)= R_close(i,:)/nanmean(R_close(i,1:4));
end

for i=1:size(R_far)
    R_close_norm(i,:)= R_close(i,:)/nanmean(R_close(i,1:4));
    % R_middle_norm(i,:)=R_middle/nanmean(R_middle(i,1:4));
    %R_close_norm(i,:)= R_close(i,:)/nanmean(R_close(i,1:4));
end

 R_far_norm_avg= nanmean(R_far_norm);
 R_far_norm_std=nanstd(R_far_norm);

 R_far_norm_lh(:,1)=R_far_norm_avg+R_far_norm_std;
 R_far_norm_lh(:,2)=R_far_norm_avg-R_far_norm_std;

 R_close_norm_avg=nanmean(R_close_norm);
 R_close_norm_std=nanstd(R_close_norm);
 R_close_norm_lh(:,1)=R_close_norm_avg+R_close_norm_std;
 R_close_norm_lh(:,2)=R_close_norm_avg-R_close_norm_std;


 R_middle_norm_avg=nanmean(R_middle_norm);
 R_middle_norm_std=nanstd(R_middle_norm);
 R_middle_norm_lh(:,1)=R_middle_norm_avg+R_middle_norm_std;
 R_middle_norm_lh(:,2)=R_middle_norm_avg-R_middle_norm_std;


x_range=[-4,-3,-2,-1,0,1,2,3,4,5,6,9,13,20,27,34,41,48,55];
f1=figure
%S=patch([x_range fliplr(x_range)],[R_far_norm_lh(:,2)' fliplr(R_far_norm_lh(:,1)')],'k','EdgeColor','k');
%hold on
plot(x_range, R_far_norm_avg,'Color','k','LineWidth',6);
hold on
%S=patch([x_range fliplr(x_range)],[R_middle_norm_lh(:,2)' fliplr(R_middle_norm_lh(:,1)')],'c','EdgeColor','c');
plot(x_range,R_middle_norm_avg,'Color','c','LineWidth',6);
%S=patch([x_range fliplr(x_range)],[R_close_norm_lh(:,2)' fliplr(R_close_norm_lh(:,1)')],'r','EdgeColor','r');
plot(x_range, R_close_norm_avg,'Color','r','LineWidth',6);
legend({'>1000μm','600~1000μm','<600μm'},'FontSize',17);
title('Ephys Far vs Close Shanks','FontSize',36);
xlabel('Days','FontSize',36)
ylabel('Normalized FR ','FontSize',36);
xticks(x_range);
xticklabels(x_range);
%xticklabels(week_info);
ax=gca;
ax.XAxis.FontSize = 22;
ax.YAxis.FontSize = 22;

save_path='Z:\xl_stroke\yifu_Digigait\SpikeSorting\Locations Analysis\Locations\Results';
saveas(f1,fullfile(save_path,'Far vs Close Shanks.fig'));
saveas(f,fullfile(save_path,'Distance Distribution.fig'));

