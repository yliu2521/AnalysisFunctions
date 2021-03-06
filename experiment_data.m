%            [1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48   49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64   65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80   81   82   83   84   85   86   87   88   89   90   91   92   93   94   95   96   97   98   99   100  101  102  103  104  105  106  107  108  109  110  111  112  113  114  115  116  117  118  119  120  121  122  123  124  125  126  127  128  129  130  131  132  133  134  135  136  137  138  139  140  141  142  143  144  145  146  147  148  149  150]
% linear track data
a = 10^(-3)*[9.05 15.4 18.7 20.3 21.3 22.1 20.7 19.9 21.1 20.9 20.5 19.2 19.2 18.4 18.1 16.2 16.3 13.1 14.3 14.3 11.3 11.3 11.9 11.5 11.2 10.1 9.59 9.74 7.94 8.62 8.48 8.72 6.87 8.87 6.73 6.68 6.14 5.90 6.05 6.34 5.22 5.27 4.88 5.90 5.51 4.67 4.93 5.08 5.08 4.25 4.01 4.25 4.10 3.96 3.91 3.96 4.20 4.15 4.06 3.86 3.91 3.48 3.48 3.48 3.82 3.18 4.01 2.41 2.85 3.09 2.51 2.41 3.09 2.21 3.09 2.56 2.99 3.19 3.53 2.46 2.65 1.92 2.22 2.41 2.31 2.95 1.78 2.36 2.27 2.37 2.80 1.98 1.34 2.56 1.98 2.42 2.71 1.64 1.83 1.74 1.59 2.22 2.03 2.32 1.64 2.18 2.03 1.88 1.45 2.61 1.01 1.69 1.25 2.96 1.74 1.94 1.79 1.45 0.96 1.16 1.26 1.40 1.74 1.40 1.26 1.50 1.50 2.04 0.97 1.60 0.87 1.31 1.45 1.50 1.45 1.36 1.89 1.60 1.26 1.60 1.60 1.55 1.26 1.17 2.19 1.21 1.07 1.02 0.78 0.97];
% open field data
b = 10^(-3)*[2.04 4.38 6.33 6.82 9.40 9.15 9.88 10.1 11.9 10.9 12.0 12.0 13.8 12.6 12.6 13.0 13.9 13.1 11.6 11.2 12.0 10.9 10.5 10.3 10.4 9.40 11.4 10.6 10.3 9.93 8.91 8.71 8.37 9.25 8.42 8.03 7.84 7.89 7.98 8.13 6.91 6.57 7.30 7.45 7.01 7.01 7.06 6.33 6.23 5.99 5.45 5.65 6.67 5.74 5.65 4.92 4.82 5.11 5.06 4.92 5.26 4.43 5.26 5.60 4.63 4.43 4.48 4.43 3.89 3.65 4.43 4.19 4.14 3.51 3.36 4.28 3.07 3.31 3.12 2.82 3.55 2.92 3.26 3.31 2.39 3.36 3.21 2.58 2.48 2.68 3.07 2.58 3.12 2.73 3.31 3.12 2.53 2.87 2.43 2.19 2.92 2.09 2.92 2.29 2.19 2.34 2.43 2.48 1.75 2.87 2.19 2.04 1.85 2.14 1.46 2.43 2.24 1.95 2.00 2.29 1.95 2.00 1.36 2.14 2.04 1.95 1.70 1.70 2.19 2.00 1.75 1.66 1.61 1.41 1.70 1.36 1.41 1.75 1.56 1.46 0.78 1.75 1.61 1.36 1.51 1.75 1.75 1.56 1.66 1.61];
c = 0.1*[1:150];
% figure
subplot(1,2,1)
% bar(0.1*(1:150),a,'histc');ylim([0 0.05])
loglog(c,a,'.')
% set(gca,'XScale','log')
xlabel('Event Step Size')
ylabel('Probability')
x = c(11:end);
v = polyfit(log10(x),log10(a(11:end)),1);
y = 10^v(2)*x.^v(1);
hold on;
loglog(x,y)
str = ['p = ',num2str(v(1))];
text(min(x),max(y),str)
title('Open field')
subplot(1,2,2)
loglog(c,b,'.')
xlabel('Event Step Size')
ylabel('Probability')
x = c(31:end);
v = polyfit(log10(x),log10(b(31:end)),1);
y = 10^v(2)*x.^v(1);
hold on;
loglog(x,y)
str = ['p = ',num2str(v(1))];
text(min(x),max(y),str)
title('Linear track')

% sum(a)