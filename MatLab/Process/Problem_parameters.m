clc
clear 
close all
%%
%Used to generate data structures for our simulations
prompt = {'Tank_name', 'Tank Height(m):','Tank volume(m^3):'...
           'Heat capacity of the fluid(J/kg°K):', 'Density of the fluid(kg/m^3):'...
           'Thermal losses coefficient(W m^{−2} K^{−1})'};
dlgtitle = 'EWH Water Tank paramaters';
dims = [1 45];
definput = {'Tank Without name', '1.12','0.112', '4182', '1000', '1.5' };
opts.Interpreter = 'tex';

answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
%Now store it in a data structure
Tank = struct('Vol', {}, 'H', {}, 'Cv', {},'Rho', {}, 'Dc', {}, 'UL', {}, 'UL_', {}, 'A', {});
Tank(1).H = str2double(answer(2)); 
Tank(1).Vol = str2double(answer(3)) ; 
Tank(1).A = Tank.Vol / Tank.H;
Tank(1).Cv = str2double(answer(4)); 
Tank(1).Rho = str2double(answer(5)); 
Tank(1).Dc = 1.8/(Tank.Rho*Tank.Cv); 
Tank(1).UL = 4*str2double(answer(6))*sqrt(pi/Tank.A)/(Tank.Rho*Tank.Cv); 
file_name = strcat(answer(1),'.mat');


%%
%Now the heating elements
prompt = {'Efficiency]0;1[:', 'Electrical power delivered by each heating element(W):'...
           'Positions of the He(m):', 'Positions of the thermos linked at each He(m):'};
dlgtitle =" EWH's Heating Elements paramaters";
dims = [1 45];
definput = {'0.95', '6000','0.2975; 0.7735', '0.2975; 0.7735'};
opts.Interpreter = 'tex';

answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
HeatElem = struct('n_eff', {}, 'Power', {}, 'Positions', {}, 'Thermos', {},'N', {});
HeatElem(1).n_eff = str2double(answer(1)); %Efficiency of the heating elements
HeatElem(1).Power = str2double(answer(2)); %Watt % Electrical power delivered by each heating element
HeatElem(1).Positions = str2num(cell2mat(answer(3)));
HeatElem(1).Thermos = str2num(cell2mat(answer(4)));
HeatElem(1).N = min([size(HeatElem.Positions, 1) size(HeatElem.Thermos, 1)]);
%%
save( char(file_name), 'Tank', 'HeatElem')
