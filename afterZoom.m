function afterZoom(~,objEvent,x,y)
% Change axes y-limits to match only what data is shown
hAx = objEvent.Axes;
yInView = y(x >= hAx.XLim(1) && x <= hAx.XLim(2));
hAx.YLim = ([min(yInView) max(yInView)]);
end 