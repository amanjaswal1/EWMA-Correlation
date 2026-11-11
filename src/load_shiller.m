function d = load_shiller(file)
%LOAD_SHILLER  Load data/shiller_monthly.csv and build monthly returns.
%
%   d = load_shiller()
%   d = load_shiller(file)
%
%   d.year     decimal year of each return (month end ~ t + 1/12)
%   d.date     'YYYY-MM' labels
%   d.stock    S&P 500 total return (price + dividend/12)
%   d.bond     10-year par-bond return from the yield (see BOND_RETURN)
%   d.X        [stock bond], T x 2

    if nargin < 1
        here = fileparts(mfilename('fullpath'));
        file = fullfile(here, '..', 'data', 'shiller_monthly.csv');
    end
    fid = fopen(file, 'r');
    if fid < 0
        error('load_shiller:file', 'Cannot open %s', file);
    end
    fgetl(fid);                                    % header
    C = textscan(fid, '%s %f %f %f', 'Delimiter', ',');
    fclose(fid);
    dates = C{1}; P = C{2}; D = C{3}; y = C{4} / 100;

    yr = cellfun(@(s) str2double(s(1:4)), dates);
    mo = cellfun(@(s) str2double(s(6:7)), dates);
    d.stock = (P(2:end) + D(2:end) / 12) ./ P(1:end-1) - 1;
    d.bond  = bond_return(y(1:end-1), y(2:end));
    d.year  = yr(2:end) + (mo(2:end) - 1) / 12;
    d.date  = cellfun(@(s) s(1:7), dates(2:end), 'UniformOutput', false);
    d.X     = [d.stock, d.bond];
    d.yield = y(2:end);
end
