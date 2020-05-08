clear all
global tab j_list

% set bootstrap to 1 to run boot strap analysis, 0 otherwise
bootstrap = 0;

if bootstrap == 1
   nruns=1000 ;
else
   nruns=1; 
end

x_out=[]
tab = readtable("mO2dat_for_model_fitting.csv");
i=1
for i = 1:nruns
if bootstrap == 1    
    j_list=     randsample(15,15,true);
else
    j_list=  1:15;
end

% x0 = (K,G,b0,b1)
x0= [1.5  , 1.7328 ,   .229 ,   0.08];
[x_out_temp val]= fminsearch(@fitO2depletion,x0)
x_out(i,:)=x_out_temp;
end

% x_out is the fitted parameters, if bootstrap analysis is run, a matrix of
% 1000 parameter sets is generated
