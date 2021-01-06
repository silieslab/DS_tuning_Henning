
function obj1 = Plot_circHist(DATA, COLOR, AXIS, rLIM, nbins)

            obj1=CircHist(convert_angle(angle([DATA])),nbins,'parent', AXIS);
            
            obj1.colorBar=COLOR;
            obj1.avgAngH.LineStyle = '--'; % make average-angle line dashed
            obj1.avgAngH.LineWidth = 1; % make average-angle line thinner
            obj1.colorAvgAng = [.5 .5 .5];
            obj1.setRLim(rLIM);
            obj1.scaleBarSide = 'right'; % draw rho-axis on the right side of the plot
            obj1.polarAxs.ThetaZeroLocation = 'right'; % rotate the plot to have 0� on the right side
            delete(obj1.rH)
%             obj1.drawArrow(obj1.avgAng, obj1.r * range(rLIM), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'k')
%             delete(obj1.scaleBar)
%             delete(obj1.thetaLabel)


            
end 