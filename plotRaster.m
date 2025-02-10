function [f]=plotRaster(spikeTrains, timeWindow, scaleBarDuration)
    % plotRaster generates a raster plot for the given spike trains
    %
    % Inputs:
    %   spikeTrains - cell array where each cell contains a vector of spike times (in seconds)
    %   timeWindow - [start_time, end_time] (in seconds), defines the time window for the plot
    %   scaleBarDuration - duration of the scale bar in seconds
    %
    % Example usage:
    %   spikeTrains = {randi([0, 1000], [50, 1]) / 1000, randi([0, 1000], [60, 1]) / 1000, randi([0, 1000], [70, 1]) / 1000};
    %   timeWindow = [0, 1]; % Time window in seconds
    %   scaleBarDuration = 0.1; % 100 ms scale bar in seconds
    %   plotRaster(spikeTrains, timeWindow, scaleBarDuration);

    if nargin < 3
        error('Three inputs required: spikeTrains, timeWindow, and scaleBarDuration');
    end

    % Plot parameters
    numNeurons = length(spikeTrains);
    yOffset = 0.3; % Offset for each spike in the raster plot

    f=figure;
    pictureSize=[200 200 2000 1200];
    set(gcf,'Position',pictureSize);
    set(gcf,'color','w')
    hold on;
    box off
    % Plot each spike train
    for neuronIdx = 1:numNeurons
        spikes = spikeTrains{neuronIdx};
        for spikeIdx = 1:length(spikes)
            spikeTime = spikes(spikeIdx);
            if spikeTime >= timeWindow(1) && spikeTime <= timeWindow(2)
                line([spikeTime, spikeTime], [neuronIdx - yOffset, neuronIdx + yOffset], 'Color', 'k','LineWidth',2);
            end
        end
    end

    % Axis labels and title
    xlabel('Time (s)');
    ylabel('Neuron');
    title('Raster Plot');
    xlim(timeWindow);
    ylim([0, numNeurons + 1]);
    set(gca, 'YTick', 1:5:numNeurons);
    axis off
    % Add scale bar
    scaleBarX = [timeWindow(1) + (timeWindow(2) - timeWindow(1)) * 0.05, ...
                 timeWindow(1) + (timeWindow(2) - timeWindow(1)) * 0.05 + scaleBarDuration];
    scaleBarY = [numNeurons + 0.5, numNeurons + 0.5];
    line(scaleBarX, scaleBarY, 'Color', 'k', 'LineWidth', 2);
    text(mean(scaleBarX), numNeurons + 1, sprintf('%.2f s', scaleBarDuration), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    hold off;
end