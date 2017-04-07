data = load('data/datarockMA4.txt');
[m,n] = size(data);
cl36AMS = data(:,n-1) ; % sample concentration in [36Cl] measured by AMS
sig_cl36AMS = data(:,n) ; % uncertainty on [36Cl] AMS measurements

avg_sig = mean(sig_cl36AMS);
%minx = cl36AMS(:, 1) - 4*sig_cl36AMS(:, 1) ;
%maxx = cl36AMS(:, 1) + 4*sig_cl36AMS(:, 1) ;

%cumminx = minx(1, 1)
%cummaxx = maxx(m, 1) + 4*sig_cl36AMS(m, 1)

minx = floor(min(cl36AMS) - 4 * avg_sig);
maxx = ceil(max(cl36AMS) + 4 * avg_sig);
x = (minx : maxx)';

cumy = zeros(maxx - minx + 1, 1);
cumavy = zeros(maxx - minx + 1, 1);
hold on;
yyaxis left;
for i = 1:m
    y = pdf('Normal', x, cl36AMS(i,1), sig_cl36AMS(i, 1));
    %y = pdf('Normal', x, cl36AMS(i,1), avg_sig);
    %plot(x, y, 'r')
    cumy = cumy + y; 
    cumavy = cumavy + pdf('Normal', x, cl36AMS(i,1), avg_sig); 
end

% plot(x, cumy, 'b');
plot(x, cumavy, 'm');

[pks, locs, w, p] = findpeaks(cumavy);
for j = 1:size(pks) 
    plot([x(locs(j)) x(locs(j))], ylim,'g-')
end
yyaxis right;
% plot data
h = data(:,n-3);
% error bar coordinates
errorX = [cl36AMS(:)-sig_cl36AMS(:) cl36AMS(:)+sig_cl36AMS(:)] ;
errorX = errorX' ;
errorY = [h(:) h(:)] ;
errorY = errorY' ;

plot(cl36AMS,h/100,'k.',errorX,errorY/100,'k-');
