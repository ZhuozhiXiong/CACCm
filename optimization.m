x0 = [0.54 0.92 1.11 0.5 1];

options4 = optimoptions('fmincon','Display','iter','PlotFcn',...
    'optimplotfvalconstr');

[solution,objectiveValue] = fmincon(@objectiveFcn,x0,[],[],[],[],...
    zeros(size(x0)),repmat(2,size(x0)),[],options4);