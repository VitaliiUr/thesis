clear all

filename="Py_30mev_scs.dat";
column_number = 3;

data02=dlmread("2.dat4.0.2.30mev","",1,0);
data12=dlmread("2.dat4.1.2.30mev","",1,0);
data22=dlmread("2.dat4.2.2.30mev","",1,0);
data32=dlmread("2.dat4.3.2.30mev","",1,0);
data42=dlmread("2.dat4.4.2.30mev","",1,0);


Pi=3.14159265359;

for i=1:length(data02)
  data02(i,1)=round(data02(i,1)*180/Pi);
endfor

data2(:,1)=data02(1:181,1);
data2(:,2)=data02(1:181,column_number);
data2(:,3)=data12(1:181,column_number);
data2(:,4)=data22(1:181,column_number);
data2(:,5)=data32(1:181,column_number);
data2(:,6)=data42(1:181,column_number);


datfile = fopen(filename,'wt');
fprintf(datfile, 'DEG\tLO\tNLO\tN2LO\tN3LO\tN4LO\n');
fclose(datfile)

dlmwrite(filename, data2,'delimiter', '\t', '-append');

clear all

filename="Py_100mev_scs.dat";
column_number = 3;

data02=dlmread("2.dat4.0.2.100mev","",1,0);
data12=dlmread("2.dat4.1.2.100mev","",1,0);
data22=dlmread("2.dat4.2.2.100mev","",1,0);
data32=dlmread("2.dat4.3.2.100mev","",1,0);
data42=dlmread("2.dat4.4.2.100mev","",1,0);


Pi=3.14159265359;

for i=1:length(data02)
  data02(i,1)=round(data02(i,1)*180/Pi);
endfor

data2(:,1)=data02(1:181,1);
data2(:,2)=data02(1:181,column_number);
data2(:,3)=data12(1:181,column_number);
data2(:,4)=data22(1:181,column_number);
data2(:,5)=data32(1:181,column_number);
data2(:,6)=data42(1:181,column_number);


datfile = fopen(filename,'wt');
fprintf(datfile, 'DEG\tLO\tNLO\tN2LO\tN3LO\tN4LO\n');
fclose(datfile)

dlmwrite(filename, data2,'delimiter', '\t', '-append');
