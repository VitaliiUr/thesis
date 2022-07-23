dataNeuNC02 = dlmread("neutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.0.2.txt","",0,0);
dataNeuNC12 = dlmread("neutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.1.2.txt","",0,0);
dataNeuNC22 = dlmread("neutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.2.2.txt","",0,0);
dataNeuNC32 = dlmread("neutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.3.2.txt","",0,0);
dataNeuNC41 = dlmread("neutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.4.1.txt","",0,0);
dataNeuNC42 = dlmread("neutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.4.2.txt","",0,0);
dataNeuNC43 = dlmread("neutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.4.3.txt","",0,0);
dataNeuNC44 = dlmread("neutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.4.4.txt","",0,0);
dataNeuNC52 = dlmread("neutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.5.2.txt","",0,0);

dataAntiNeuNC02 = dlmread("antineutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.0.2.txt","",0,0);
dataAntiNeuNC12 = dlmread("antineutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.1.2.txt","",0,0);
dataAntiNeuNC22 = dlmread("antineutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.2.2.txt","",0,0);
dataAntiNeuNC32 = dlmread("antineutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.3.2.txt","",0,0);
dataAntiNeuNC41 = dlmread("antineutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.4.1.txt","",0,0);
dataAntiNeuNC42 = dlmread("antineutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.4.2.txt","",0,0);
dataAntiNeuNC43 = dlmread("antineutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.4.3.txt","",0,0);
dataAntiNeuNC44 = dlmread("antineutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.4.4.txt","",0,0);
dataAntiNeuNC52 = dlmread("antineutrino_NC_deuteron_xs_results_2019_from_table_REGIONS.SMS.5.2.txt","",0,0);

dataAntiNeuCC02 = dlmread("antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.0.2.txt","",0,0);
dataAntiNeuCC12 = dlmread("antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.1.2.txt","",0,0);
dataAntiNeuCC22 = dlmread("antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.2.2.txt","",0,0);
dataAntiNeuCC32 = dlmread("antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.3.2.txt","",0,0);
dataAntiNeuCC41 = dlmread("antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.1.txt","",0,0);
dataAntiNeuCC42 = dlmread("antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.2.txt","",0,0);
dataAntiNeuCC43 = dlmread("antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.3.txt","",0,0);
dataAntiNeuCC44 = dlmread("antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.4.4.txt","",0,0);
dataAntiNeuCC52 = dlmread("antineutrino_CC_deuteron_xs_results_2019_from_table.SMS.5.2.txt","",0,0);

# E = 20
num=11;
display(["E = ", num2str(dataNeuNC42(num,1))])

display("N4LO-N4LO+ spread at E = 20 MeV (%)")
display("neutrino NC")
abs(dataNeuNC42(num,3)-dataNeuNC52(num,3))/(dataNeuNC42(num,3)+dataNeuNC52(num,3))/0.5*100

display("antineutrino NC")
abs(dataAntiNeuNC42(num,3)-dataAntiNeuNC52(num,3))/(dataAntiNeuNC42(num,3)+dataAntiNeuNC52(num,3))/0.5*100

display("antineutrino CC")
abs(dataAntiNeuCC42(num,3)-dataAntiNeuCC52(num,3))/(dataAntiNeuCC42(num,3)+dataAntiNeuCC52(num,3))/0.5*100
""

display("cut-off spread at E = 20 MeV (%)")
display("neutrino NC")
abs(max([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)])-...
	min([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))/...
	abs(mean([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))*100

display("antineutrino NC")
abs(max([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)])-...
	min([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))/...
	abs(mean([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))*100

display("antineutrino CC")
abs(max([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)])-...
	min([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))/...
	abs(mean([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))*100

# E = 30
num=16;
display(["E = ", num2str(dataNeuNC42(num,1))])
display("N4LO-N4LO+ spread at E = 30 MeV (%)")
display("neutrino NC")
abs(dataNeuNC42(num,3)-dataNeuNC52(num,3))/(dataNeuNC42(num,3)+dataNeuNC52(num,3))/0.5*100

display("antineutrino NC")
abs(dataAntiNeuNC42(num,3)-dataAntiNeuNC52(num,3))/(dataAntiNeuNC42(num,3)+dataAntiNeuNC52(num,3))/0.5*100

display("antineutrino CC")
abs(dataAntiNeuCC42(num,3)-dataAntiNeuCC52(num,3))/(dataAntiNeuCC42(num,3)+dataAntiNeuCC52(num,3))/0.5*100
""

display("cut-off spread at E = 30 MeV (%)")
display("neutrino NC")
abs(max([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)])-...
	min([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))/...
	abs(mean([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))*100

display("antineutrino NC")
abs(max([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)])-...
	min([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))/...
	abs(mean([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))*100

display("antineutrino CC")
abs(max([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)])-...
	min([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))/...
	abs(mean([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))*100

# E = 100
num=51;
display(["E = ", num2str(dataNeuNC42(num,1))])

display("N4LO-N4LO+ spread at E = 100 MeV (%)")
display("neutrino NC")
abs(dataNeuNC42(num,3)-dataNeuNC52(num,3))/(dataNeuNC42(num,3)+dataNeuNC52(num,3))/0.5*100

display("antineutrino NC")
abs(dataAntiNeuNC42(num,3)-dataAntiNeuNC52(num,3))/(dataAntiNeuNC42(num,3)+dataAntiNeuNC52(num,3))/0.5*100

display("antineutrino CC")
abs(dataAntiNeuCC42(num,3)-dataAntiNeuCC52(num,3))/(dataAntiNeuCC42(num,3)+dataAntiNeuCC52(num,3))/0.5*100
""

display("cut-off spread at E = 100 MeV (%)")
display("neutrino NC")
abs(max([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)])-...
	min([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))/...
	abs(mean([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))*100

display("antineutrino NC")
abs(max([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)])-...
	min([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))/...
	abs(mean([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))*100

display("antineutrino CC")
abs(max([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)])-...
	min([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))/...
	abs(mean([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))*100

# E = 180
num=83;
display(["E = ", num2str(dataNeuNC42(num,1))])

display("N4LO-N4LO+ spread at E = 180 MeV (%)")
display("neutrino NC")
abs(dataNeuNC42(num,3)-dataNeuNC52(num,3))/(dataNeuNC42(num,3)+dataNeuNC52(num,3))/0.5*100

display("antineutrino NC")
abs(dataAntiNeuNC42(num,3)-dataAntiNeuNC52(num,3))/(dataAntiNeuNC42(num,3)+dataAntiNeuNC52(num,3))/0.5*100

display("antineutrino CC")
abs(dataAntiNeuCC42(num,3)-dataAntiNeuCC52(num,3))/(dataAntiNeuCC42(num,3)+dataAntiNeuCC52(num,3))/0.5*100
""

display("cut-off spread at E = 180 MeV (%)")
display("neutrino NC")
abs(max([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)])-...
	min([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))/...
	abs(mean([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))*100

display("antineutrino NC")
abs(max([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)])-...
	min([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))/...
	abs(mean([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))*100

display("antineutrino CC")
abs(max([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)])-...
	min([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))/...
	abs(mean([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))*100

NeuNCchir = [];
AntiNeuNCchir = [];
AntiNeuCCchir = [];

NeuNCcut = [];
AntiNeuNCcut = [];
AntiNeuCCcut = [];
energy = [];
for i = 1:85
	% if(i!=2 && i!= 4)
		num = i;
		energy = [energy;dataAntiNeuNC42(i,1)];
		NeuNCchir = [NeuNCchir; abs(dataNeuNC42(num,3)-dataNeuNC52(num,3))/(dataNeuNC42(num,3)+dataNeuNC52(num,3))/0.5*100];

		AntiNeuNCchir = [AntiNeuNCchir; abs(dataAntiNeuNC42(num,3)-dataAntiNeuNC52(num,3))/(dataAntiNeuNC42(num,3)+dataAntiNeuNC52(num,3))/0.5*100];

		AntiNeuCCchir = [AntiNeuCCchir; abs(dataAntiNeuCC42(num,3)-dataAntiNeuCC52(num,3))/(dataAntiNeuCC42(num,3)+dataAntiNeuCC52(num,3))/0.5*100];

		% NeuNCcut = [NeuNCcut; abs(max([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)])-...
		% 	min([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))/...
		% 	abs(mean([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))*100];

		% AntiNeuNCcut = [AntiNeuNCcut; abs(max([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)])-...
		% 	min([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))/...
		% 	abs(mean([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))*100];

		% AntiNeuCCcut = [AntiNeuCCcut; abs(max([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)])-...
		% 	min([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))/...
		% 	abs(mean([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))*100];
		NeuNCcut = [NeuNCcut; abs(max([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)])-...
			min([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))/...
			abs(max([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)])+...
			min([dataNeuNC41(num,3),dataNeuNC42(num,3),dataNeuNC43(num,3),dataNeuNC44(num,3)]))/0.5*100];

		AntiNeuNCcut = [AntiNeuNCcut; abs(max([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)])-...
			min([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))/...
			abs(max([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)])+...
			min([dataAntiNeuNC41(num,3),dataAntiNeuNC42(num,3),dataAntiNeuNC43(num,3),dataAntiNeuNC44(num,3)]))/0.5*100];

		AntiNeuCCcut = [AntiNeuCCcut; abs(max([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)])-...
			min([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))/...
			abs(max([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)])+...
			min([dataAntiNeuCC41(num,3),dataAntiNeuCC42(num,3),dataAntiNeuCC43(num,3),dataAntiNeuCC44(num,3)]))/0.5*100];
	% endif
endfor

plot(energy,[NeuNCchir,AntiNeuNCchir,AntiNeuCCchir])
figure(2)
plot(energy,[NeuNCcut,AntiNeuNCcut,AntiNeuCCcut])

data = [energy,NeuNCchir,AntiNeuNCchir,AntiNeuCCchir,NeuNCcut,AntiNeuNCcut,AntiNeuCCcut];

datfile = fopen("errors.dat",'wt');
fprintf(datfile, 'E\tnuNC_chiral\tnubarNC_chiral\tnubarCC_chiral\tnuNC_cutoff\tnubarNC_cutoff\tnubarCC_cutoff\n');
fclose(datfile)

dlmwrite("errors.dat", data,'delimiter', '\t', '-append', 'precision', 10);
