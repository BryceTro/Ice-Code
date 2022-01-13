%LHS 1140 b Comparison Plots

figure
subplot(2,3,1)
T1 = readtable('ResearchCandidateList.xlsx');
vars = T1.Properties.VariableNames;
plot(T1{:,5},'-o')
Ax = gca;
Ax.XTick = 1:numel(T1{:,1});
Ax.XTickLabel = T1{:,1};
xlabel(vars{1})
ylabel(vars{5})
xtickangle(45)
grid on;

subplot(2,3,2)
T2 = readtable('ResearchCandidateList.xlsx');
vars = T2.Properties.VariableNames;
plot(T2{:,6},'-o')
Ax = gca;
Ax.XTick = 1:numel(T2{:,1});
Ax.XTickLabel = T2{:,1};
xlabel(vars{1})
ylabel(vars{6})
xtickangle(45)
grid on;

subplot(2,3,3)
T3 = readtable('ResearchCandidateList.xlsx');
vars = T3.Properties.VariableNames;
plot(T3{:,7},'-o')
Ax = gca;
Ax.XTick = 1:numel(T3{:,1});
Ax.XTickLabel = T3{:,1};
xlabel(vars{1})
ylabel(vars{7})
xtickangle(45)
grid on;

subplot(2,3,4)
T4 = readtable('ResearchCandidateList.xlsx');
vars = T4.Properties.VariableNames;
plot(T4{:,8},'-o')
Ax = gca;
Ax.XTick = 1:numel(T1{:,1});
Ax.XTickLabel = T1{:,1};
xlabel(vars{1})
ylabel(vars{8})
xtickangle(45)
grid on;

subplot(2,3,5)
T5 = readtable('ResearchCandidateList.xlsx');
vars = T5.Properties.VariableNames;
plot(T5{:,9},'-o')
Ax = gca;
Ax.XTick = 1:numel(T5{:,1});
Ax.XTickLabel = T5{:,1};
xlabel(vars{1})
ylabel(vars{9})
xtickangle(45)
grid on;

subplot(2,3,6)
T6 = readtable('ResearchCandidateList.xlsx');
vars = T6.Properties.VariableNames;
plot(T6{:,10},'-o')
Ax = gca;
Ax.XTick = 1:numel(T6{:,1});
Ax.XTickLabel = T6{:,1};
xlabel(vars{1})
ylabel(vars{10})
xtickangle(45)
grid on;


