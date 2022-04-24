---
title: "Power Analysis"
author: "Chun-Jie Liu"
date: "2022-04-22"
---

Power analysis is an important aspect of experimental design. It allows us to determine the sample size required to detect an effect of a given size with a given degree of confidence. Conversely, it allows us to determine the probability of detecting an effect of a given size with a given size with a given level of confidence, under sample size contraints.

1. sample size
2. effect size
3. significance level = P(Type I error) = probability of of finding an effect that is not there
4. power = 1 - P(Type II error) = probability of finding an effect that is there

Given any three, we can determine the forth.

The significance level defaults to 0.05. Therefore, to calculate the significance level, given an effect size, sample size, and power, use the option "sig.level = NULL".