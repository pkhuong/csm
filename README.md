Safely-rounded Confidence Sequence Method
=========================================

An implementation of [Ding, Gandy, and Hahn](https://arxiv.org/abs/1611.01675)'s
Confidence Sequence Method for dynamic termination of Binomial tests on
Monte Carlo simulations, as well as code to compute Bayesian credible
intervals on the underlying success rate.

The code in C, Common Lisp, and Python should give bitwise identical
results, and the results should always be safely (conservatively)
rounded; the implementation of CSM in floating point should not cause
any false positive. N.B., The reason I care about round-off errors
isn't that I expect them to make a big difference, but simply that I
want peace of mind, and not have to wonder if weird numerical issues
are causing whatever my system reports.

See [this blog post](https://www.pvk.ca/Blog/2018/07/06/testing-slo-type-properties-with-the-confidence-sequence-method/) for more information.

I'd particularly like to receive suggestions or pull requests to
improve the usability of the tools. I'm thinking improvements to the
report format, graphical reports, and CLIs tools that are easier to
insert in testing pipelines or in ad hoc shell scripts.
