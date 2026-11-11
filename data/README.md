# Data

`shiller_monthly.csv`: monthly U.S. stock and bond data, April 1953 to June 2023 (843 rows).

| Column | Meaning |
|---|---|
| `date` | First day of the month |
| `sp500_price` | S&P 500 composite index, **monthly average** of daily closes |
| `sp500_dividend_annual` | Dividends per share, annualised |
| `long_rate_pct` | 10-year U.S. Treasury constant-maturity yield, % per year (FRED GS10, monthly average) |

**Source:** Robert Shiller's U.S. stock market data (<http://www.econ.yale.edu/~shiller/data.htm>), as packaged by the Open Knowledge Foundation's `datasets/s-and-p-500` repository on GitHub. **License:** Open Data Commons Public Domain Dedication and License (PDDL) v1.0.

**Why this window**
- It starts in April 1953, when the 10-year Treasury series (GS10) begins. Before that, Shiller's long rate is interpolated from annual data, so month-to-month bond returns computed from it would be artificial.
- It ends in June 2023, the last month with dividend data in the packaged file.

**Refresh:** `data/fetch_shiller.m` downloads the latest packaged file and rewrites this CSV with the same filters.

**How returns are built** (`src/load_shiller.m`, `src/bond_return.m`)
- **Stock:** monthly total return $r_s = (P_t + D_t/12)/P_{t-1} - 1$.
- **Bond:** the return on a 10-year par bond bought last month at yield $y_{t-1}$ and repriced at $y_t$: coupon income $y_{t-1}/12$ plus the price change. This is a standard approximation that ignores one month of roll-down.
