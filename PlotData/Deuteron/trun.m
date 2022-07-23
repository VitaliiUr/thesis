data30=dlmread("data/Trunc_VU30MeV_cross_sms.dat","",1,0);
% data30=dlmread("data/Thruns_VU30MeV.dat","",1,0);
data30r = data30(1:180,1);
for i = 1:180
	data30r(i,2) = abs(data30(i,6)-data30(361-i,6));
endfor

[el30,num30] = max(data30r(:,2));

data30r(num30,2)/(mean([data30(num30,6),data30(361-num30,6)]))*100;

disp(cstrcat("For 30 MeV max error is ",num2str(ans)," %"))

data100=dlmread("data/Trunc_VU100MeV_cross_sms.dat","",1,0);
% data100=dlmread("data/Thruns_VU100MeV.dat","",1,0);
data100r = data100(1:180,1);
for i = 1:180
	% data100r(i,2:6) = abs(data100(i,2:6)-data100(361-i,2:6));
	data100r(i,2) = abs(data100(i,6)-data100(361-i,6));
endfor

% [el100,num100] = max(data100r(:,6));
% data100r(num100,6)/(mean([data100(num100,6),data100(361-num100,6)]))*100;
[el100,num100] = max(data100r(:,2));
data100r(num100,2)/(mean([data100(num100,6),data100(361-num100,6)]))*100;
disp(cstrcat("For 100 MeV max error is ",num2str(ans)," %"))