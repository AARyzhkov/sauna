---
icon: question
---

# Uncertainty

The uncertainty is a numerical measure how well a quantity is known and defined as "an interval having a stated level of confidence" according to the ISO 1993 terms.&#x20;

It shall be noted that the uncertainty term shall not be confused with the error term. The error term implies the difference between the so-called "true" value and the measurement while uncertainty relates to statistics having probability-related origin. At the same time, the "true" value cannot be known, only with a certain level of accuracy, therefore a reference value is chosen, usually a well-defined experiment or an average.&#x20;

To make the difference clearer, here is an example: one has a critical experiment with a multiplication factor of unity (assume it is the "true" value); the result of a stochastic simulation is 0.99980 with a stochastic uncertainty of 30 pcm, i.e. 0.99980+/- 0.00030. The error here is \~20 pcm, but the result contains an uncertainty, hence the error is 20 pcm with a statistical uncertainty of the error of 30 pcm, i.e. 20 +/- 30. That is, the error is a difference between two values while the uncertainty is complaint with the rules of statistics with the corresponding confidence level.&#x20;

The uncertainty is usually defined in a number of standard deviations $$\sigma$$, and, in nuclear data field, it is widely accepted in terms of 1-$$\sigma$$, which is equal to the 68% confidence interval (see the 68-95-99.7 rule). In this case, the uncertainty is equal to the standard deviation or the square root of the dispersion, which is aligned with that covariance matrices use dispersions (on the diagonal if symmetric) and covariances by its definition, which are considered in the following subsections.

