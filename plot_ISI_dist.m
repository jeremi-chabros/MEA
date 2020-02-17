%% get ISIs; plot CDF function and theoretical poission CDF function based
% on the mean ISI

% mu = 22;
% sigma = 1;
% pd = makedist('Normal','mu',mu,'sigma',sigma);
% % Define the input vector x to contain the values at which to calculate the cdf.
% 
% x = linspace(22-5,22+5,100);
% % Compute the cdf values for the standard normal distribution at the values in x.
% 
% y = cdf(pd,x);
% figure;
% plot(y)
% 
% %% Create a Poisson distribution object with the rate parameter, ?, equal to 2.
% % 
% % figure; plot(y)
% % figure; plot(pdf(pd,x))
% % 
% figure;
% ISIs = normrnd(10,2.5,100,1);
% histogram(ISIs,'Normalization','probability')
% aesthetics
% hold on
% 
% lambda = mean(ISIs);
% pd = makedist('Poisson','lambda',lambda);
% % Define the input vector x to contain the values at which to calculate the cdf.
% maxISI = 30;
% minISI = 1;
% x = [minISI:1:maxISI];
% % xlim([0 10])
% % Compute the cdf values for the Poisson distribution at the values in x.
% 
% y = cdf(pd,x);
% 
% plot(pdf(pd,x))
% 
% mu = mean(ISIs);
% sigma = std(ISIs);
% pd = makedist('Normal','mu',mu,'sigma',sigma);
% plot(pdf(pd,x))
% % lambda = [1:1:20];
% % y = 1 - poisscdf(2,lambda);
% % pd = makedist('Normal','mu',mu,'sigma',sigma);
% % % Define the input vector x to contain the values at which to calculate the cdf.
% % 
% % x = linspace(22-5,22+5,100);
% % % Compute the cdf values for the standard normal distribution at the values in x.
% % 
% % y = poisscdf(pd,x);
% % figure;
% % plot(y)
% % 

% load data set parameters
load('200127_FTDOrg_GrpD_5B_Slice2_mSpikes_3.mat')
if ~exist('fs')
    fs = 25000;
end

spiketimes = find(full(mSpikes(:,1))==1);
ISIs = ( spiketimes(2:end,1) - spiketimes(1:end-1,1) )  / (fs/1000);
figure
h = histogram(sqrt(ISIs),'Normalization','probability')

% peak = 1 + h.BinWidth * x(find(h.Values == max(h.Values)));

aesthetics
xticklabels(xticks.^2)
xtickformat(    )
xlabel('Inter-spike interval (ms)')
ylabel('{p       }')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'FontName','Arial','FontSize',14)
% ax = gca;
% ax.FontSize = 12
% histogram(ISIs)
mean(sqrt(ISIs))

lambda = mean(sqrt(ISIs));
% lambda = peak;
pd = makedist('Poisson','lambda',lambda);
% Define the input vector x to contain the values at which to calculate the cdf.
% maxISI = max(sqrt(ISIs));
% minISI = min(sqrt(ISIs));
x = [1:1:max(xticks)];
% xlim([0 10])
% Compute the cdf values for the Poisson distribution at the values in x.

%y = cdf(pd,x);
hold on

mu =  mean(sqrt(ISIs));
% mu = peak;
 sigma = std(sqrt(ISIs));

pdn = makedist('Normal','mu',mu,'sigma',sigma);
plot(pdf(pdn,x)/2,'LineWidth',2,'Color','b','LineStyle',':')

plot(pdf(pd,x)/2,'LineWidth',2,'Color','r','LineStyle',':')


%cdf plot 

figure
plot(cdf(pdn,x),'LineWidth',2,'Color','b','LineStyle',':')
hold on
plot(cdf(pd,x),'LineWidth',2,'Color','r','LineStyle',':')
aesthetics
xticklabels(xticks.^2)
xtickformat(    )
xlabel('Inter-spike interval (ms)')
ylabel('Cumulative p')
set(get(gca,'ylabel'),'rotation',90)
set(gca,'FontName','Arial','FontSize',14)

c =  cdfplot(sqrt(ISIs));  
grid off
c.LineWidth = 2;
c.Color = 'k';
c.Parent.XLabel.String = 'Inter-spike interval (ms)';
legend('Poisson CDF','Normal CDF','Empirical CDF','Location','best');
% sq_sorted = sort(sqrt(ISIs));
% cdf(sq_sorted)
