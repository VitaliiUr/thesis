clear all

filename="CROSS_SingCur_140mev.dat";
column_number = 3;

data42=dlmread("data/0.dat1.4.2.140mev","",2,0);
data422=dlmread("data/2.dat1.4.2.140mev","",2,0);
%dataAV18=dlmread("data/2.dat1.av18.100mev","",2,0);


% Pi=3.14159265359;

for i=1:length(data42)
  data42(i,1)=round(data42(i,1)*180/pi);
endfor

data2(:,1)=data42(1:181,1);
data2(:,2)=data42(1:181,column_number);
data2(:,3)=data422(1:181,column_number);
%data2(:,11)=dataAV18(1:180,column_number);


datfile = fopen(filename,'wt');
fprintf(datfile, 'DEG\tN4LO4_SNC\tN4LO_Sieg\n');
fclose(datfile)

dlmwrite(filename, data2,'delimiter', '\t', '-append');
