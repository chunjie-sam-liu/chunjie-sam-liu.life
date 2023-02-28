---
title: "Arrow, DuckDB, data.table and vroom"
author: "Chun-Jie Liu"
date: "2022-08-03"
---


## R packages for loading and wrangling tabular data

Fatest way to read and tidy tabular data then import it into an embedding database.

`vroom` is really fast to read the big tabular data (> 10 million observations) into R working environment, it is indeed best choice to load data into memory. However, it's very slow and taking very large memory to wranggle rows or columns into tidy data frame. Even the simplest query or filtering of one observation takes unexpected time.

`Apache arrow`

`DuckDB vs. SQLite`

`data.table` and `dtplyr`


## Best choice for practical demands

My problem is to load


