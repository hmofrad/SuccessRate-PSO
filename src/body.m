clear
clc
run = 50;
% Define dimension: D = 30 or 300 in Imp_GSO_Func
% Define Group size: Ps = 48 in Imp_GSO_Func
res = [];
f = 1; % 1 to 13 as benchmark function index
for i=1:run
    out = Imp_GSO_Func(f);
    res = [res;out(end)];
end
[mean(res) std(res)]

    