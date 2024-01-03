function  plot_modal_shape_animation(phi,varargin)
    if nargin == 0
        clc;clear;close all;
        disp('Running tests...');
        phi = [1+0.1*i,-1]';
        plot_modal_shape_animation(phi)
        disp('Tests completed.');
        return;
    end

    p = inputParser;
    addParameter(p, 'showtext', true, @islogical)
    addParameter(p, 'showplot', true, @islogical)
    addParameter(p, 'loc', 1:length(phi), @isnumeric)
    addParameter(p, 'f', 1, @isnumeric)
    addParameter(p, 'n', 100, @isnumeric)
    addParameter(p, 'm', 3, @isnumeric)
    parse(p, varargin{:});
    showtext = p.Results.showtext;
    showplot = p.Results.showplot;
    loc = p.Results.loc;
    f = p.Results.f;
    n = p.Results.n;
    m = p.Results.m;

    % for k1 = 1:n*m 
    %     u(k1,:) = real(phi*exp(1i*2*pi*f*k1/n));
    % end

    if showplot
        disp("The totol animation time is around "+num2str(m/f)+" seconds.")
        disp("There are "+num2str(m)+" cycles in total.")
        disp("The sampling frequency is "+num2str(f)+" Hz.")
        disp("There are "+num2str(n)+" samples in each cycle.")
    end
    % if showplot
    %     figure
    %     % animate
    %     for k1 = 1:n*m 
    %         plot(loc,u(k1,:),'-o','MarkerSize',10,'MarkerFaceColor','r')
    %         ylim([-1.5,1.5])
    %         xlim([0,3])
    %         grid on
    %         if showtext
    %             text(loc(1),u(k1,1),['\leftarrow',num2str(u(k1,1))])
    %             text(loc(2),u(k1,2),['\leftarrow',num2str(u(k1,2))])
    %         end
    %         pause(1/(n*f))
    %         if k1 ~= n*m 
    %             cla
    %         end
    %     end
    % end

    if showplot
        fig = figure;
        ax = axes('Parent', fig);
        btn = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Start Animation',...
                        'Position', [20 20 100 40], 'Callback', @startAnimation);
    end

    function startAnimation(src, event)
        for k1 = 1:n*m
            u(k1, :) = real(phi * exp(1i * 2 * pi * f * k1 / n));
            plot(ax, loc, u(k1, :), '-o', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
            ylim(ax, [-1.5, 1.5])
            % xlim(ax, [0, 3])
            grid(ax, 'on')
            if showtext
                text(ax, loc(1), u(k1, 1), ['\leftarrow', num2str(u(k1, 1))])
                text(ax, loc(2), u(k1, 2), ['\leftarrow', num2str(u(k1, 2))])
            end
            pause(1/(n*f))
            if k1 ~= n*m
                cla(ax)
            end
        end
    end

end
