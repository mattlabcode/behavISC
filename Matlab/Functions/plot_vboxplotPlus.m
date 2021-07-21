function xvals = plot_vboxplotPlus(x,w,c,s, plot_spread)
 
% xvals = plot_vboxplotPlus(x,w,c,s)
%
% x - data vector
% w - [min max], defines the width and also the position on the x-axis
% c - color vectors; [0 0 0; 1 0 0; 1 0 0] --> line, center, data points
% s - side of boxplot --> 0 - left, 1 - right
%
% xvals - x values for individual data points
%
% Description: The script plots a vertical boxplot on one half and individual 
% data points on the other half. Use set(gca,'XTick',..., 'XTickLabel', ...)
% to fix x-axis labels.
%
% ---------
%
%    Copyright (C) 2018
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% -----------------------------------------------------------------------
% B. Herrmann, Email: herrmann.b@gmail.com, 2018-02-25
 
if nargin < 2, fprintf('Error: Wrong number of inputs!\n'), return; end
if ~isvector(x), fprintf('Error: x needs to be a vector!\n'), return; end
 
x  = x(:);                        % column vector
d  = prctile(x,[0 25 50 75 100]); % percentiles
lw = 2;                           % line width
 
% get boundaries for plotting
if s == 0
    whalf   = [0 -diff(w)/2];
    whalfT  = [diff(w)/4 -diff(w)/2];
    wspread = [w(1)+diff(w)/2 w(2)];
else
    whalf   = [diff(w)/2 0];
    whalfT  = [diff(w)/2 -diff(w)/4];
    wspread = [w(1) w(2)-diff(w)/2];
end
 
% plot vertical line
line([w(1)+(w(2)-w(1))/2 w(1)+(w(2)-w(1))/2],[d(1) d(5)],'LineStyle','-','Color',c(1,:),'LineWidth',lw)
 
% plot ticks
line([w(1) w(2)]+whalfT,[d(1) d(1)],'LineStyle','-','Color',c(1,:),'LineWidth',lw)
line([w(1) w(2)]+whalfT,[d(5) d(5)],'LineStyle','-','Color',c(1,:),'LineWidth',lw)
 
% plot center part
if d(4)-d(2) ~= 0
	rectangle('Position',[w(1)+whalf(1) d(2) (w(2)+whalf(2))-(w(1)+whalf(1)) d(4)-d(2)],'FaceColor',c(2,:),'EdgeColor',c(1,:),'LineWidth',lw)
end
line([w(1) w(2)] + whalf,[d(3) d(3)],'LineStyle','-','Color',c(1,:),'LineWidth',lw)

% plot spread
if ~exist('plot_spread', 'var') || plot_spread == 1
    sp_handle = plotSpread(x,'xValues',mean(wspread),'distributionColors',{c(3,:)},'spreadWidth',diff(wspread),'distributionMarkers','o');
    set(sp_handle{1},'LineWidth',0.25,'MarkerFaceColor',c(3,:),'MarkerEdgeColor',[0 0 0],'MarkerSize',8);
    xvals = get(sp_handle{1},'XData');
elseif plot_spread == 0
    xvals = [];
end
