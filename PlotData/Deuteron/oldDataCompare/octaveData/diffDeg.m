data30 = dlmread("CROSS_30mev.dat","",1,0);
data100 = dlmread("CROSS_100mev.dat","",1,0);
data30_old = dlmread("CROSS_old_30mev.dat","",1,0);
data100_old = dlmread("CROSS_old_100mev.dat","",1,0);

d30 = [30, 60, abs(diff([data30(61,2:5),data30(61,7),data30(61,10)]))];
d100 = [100, 15, abs(diff([data100(16,2:5),data100(16,7),data100(16,10)]));\
      100, 150, abs(diff([data100(151,2:5),data100(151,7),data100(151,10)]))]

d30old = [30, 60,abs(diff(data30_old(61,2:6)))]
d100old = [100,15,abs(diff(data100_old(16,2:6)));100,150,abs(diff(data100_old(151,2:6)))]

% datfile = fopen("diffDegSMS.dat",'wt');
% fprintf(datfile, 'En\tTheta\tLO\tNLO\tN2LO\tN3LO\tN4LO\n');
% fclose(datfile)


dlmwrite("diffDegSMS.dat", [d30;d100]','delimiter', '\t', 'precision', 4);

% datfile = fopen("diffDegCMS.dat",'wt');
% fprintf(datfile, 'En\tTheta\tLO\tNLO\tN2LO\tN3LO\n');
% fclose(datfile)


dlmwrite("diffDegCMS.dat", [d30old;d100old]','delimiter', '\t', 'precision', 4);
