function h = shade_error_region(x, y, err, color)
% SHADE_ERROR_REGION  draw a semi-transparent error region around a line
%   h = shade_error_region(x, y, err, color)
%   x, y: vectors of same length
%   err: vector of same length containing error magnitude (symmetric)
%   color: RGB triple for the fill

x = x(:)'; y = y(:)'; err = err(:)';
% Lightweight diagnostics to assist debugging visibility issues
try
    any_nonzero = any(err(:) > 0);
    fprintf('shade_error_region: len=%d, err_min=%g, err_max=%g, any_nonzero=%d\n', numel(err), min(err(:)), max(err(:)), any_nonzero);
catch
    % ignore diagnostics errors
end

% Create x coordinates for the error region (forward and back)
x_polygon = [x(:)', fliplr(x(:)')];

% Create y coordinates for the error region (upper and lower bounds)
y_upper = (y(:)') + (err(:)');
y_lower = (y(:)') - (err(:)');
y_polygon = [y_upper, fliplr(y_lower)];

h = fill(x_polygon, y_polygon, color, 'LineStyle', 'none', 'Parent', gca);
set(h, 'EdgeColor', 'none');
% Slightly stronger transparency so the band is visible on many displays
set(h, 'FaceAlpha', 0.45); % Transparency level

% Prevent the filled polygon from appearing in legends
try
    set(h, 'HandleVisibility', 'off');
    if isprop(h, 'Annotation')
        li = get(h, 'Annotation');
        if isstruct(li) || ~isempty(li)
            set(get(h, 'Annotation').LegendInformation, 'IconDisplayStyle', 'off');
        end
    end
catch
    % Ignore if properties aren't available in older MATLAB versions
end
% Ensure the patch is underneath plotted lines so shading doesn't obscure them
try
    uistack(h, 'bottom');
catch
    % uistack may not be available in some contexts; ignore silently
end
try
    drawnow;
catch
    % ignore
end
end
