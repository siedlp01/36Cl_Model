
zero_sliprate = [
0 0;
5 0.073643411;
12 0;
13 0.125;
%14 0.036752137;
15 0.204081633;
16 0.037378641;
18 0;
];

sch_sliprate = [
0 0;
5 0.126666667;
12 0;
13 0.1875;
%14 0.058503401;
15 0.253164557;
16 0.093902439;
18 0;
];




hold on;

plot(zero_sliprate(:,1), zero_sliprate(:,2), '-r', 'LineWidth', 2);
plot(sch_sliprate(:,1), sch_sliprate(:,2), '-b', 'LineWidth', 2);

title('Slip-rate variation along fault')
ylabel('Slip-rate (cm/yr)');
xlabel('Velino-Magnola Fault (sites)');
xticks([5,13,14,15,16]);
xticklabels({'VE', 'MA1', 'MA2', 'MA3', 'MA4'});
hold off;


