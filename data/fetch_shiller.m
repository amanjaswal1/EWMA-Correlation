% FETCH_SHILLER  Download the packaged Shiller data and rebuild shiller_monthly.csv.
%
%   run('data/fetch_shiller.m')
%   Needs an internet connection. Keeps rows from April 1953 (start of the
%   10-year Treasury series) where dividends and the long rate are present.

url = 'https://raw.githubusercontent.com/datasets/s-and-p-500/main/data/data.csv';
here = fileparts(mfilename('fullpath'));
try
    txt = webread(url, weboptions('ContentType', 'text'));   % force raw text (not a table)
catch
    txt = urlread(url);
end
txt = char(txt);
lines = strsplit(strtrim(txt), sprintf('\n'));
fid = fopen(fullfile(here, 'shiller_monthly.csv'), 'w');
fprintf(fid, 'date,sp500_price,sp500_dividend_annual,long_rate_pct\n');
n = 0;
for i = 2:numel(lines)
    f = strsplit(strtrim(lines{i}), ',');
    % columns: Date,SP500,Dividend,Earnings,CPI,Long Interest Rate,...
    d = f{1};
    p = str2double(f{2}); div = str2double(f{3}); y = str2double(f{6});
    if strcmp(d, '1953-04-01') || (n > 0)
        if p > 0 && div > 0 && y > 0
            fprintf(fid, '%s,%.10g,%.10g,%.10g\n', d, p, div, y);
            n = n + 1;
        elseif n > 0
            break                      % stop at the first incomplete month
        end
    end
end
fclose(fid);
fprintf('wrote %d rows to data/shiller_monthly.csv\n', n);
