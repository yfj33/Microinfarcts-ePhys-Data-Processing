clc;clear;close all
animal_num=5;
Stroke_Size=[152; 165;163.5;150;250]
D2Rec=[14 10 10 0;14 nan nan nan;0 7 0 0;7 0 0 5;nan 21 10 6;];
WorstD=[4 5 3 nan;3 3 nan nan;nan 4 nan 2;3 3 6 3;3 10 3 2];
LowP=[0.45 0.58 0.358 nan;0.38 0.01 nan nan;nan 0.85 nan 0.82;0.63 0.91 0.79 0.785;0.08 0.47 0 0.42];

f1=figure
for i=1:animal_num
    
    scatter(Stroke_Size(i,1),D2Rec(i,:),'red','filled');
    hold on

    
end
xlabel('Radius(um)');
ylabel('Days');
title('Days to Recover to 80%');
f2=figure
   for i=1:animal_num
    
    scatter(Stroke_Size(i,1),WorstD(i,:),'green','filled');
    hold on

    
   end
   xlabel('Radius(um)');
ylabel('Days');
title('Worst Days');

    f3=figure
   for i=1:animal_num
    
    scatter(Stroke_Size(i,1),LowP(i,:),'blue','filled');
    hold on

    
   end
 xlabel('Radius(um)');
ylabel('Normalized FR');
title('Lowest FR');

