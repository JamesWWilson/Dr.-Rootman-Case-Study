# Dr.-Rootman-Eye-Polygon-Case-Study
UCLA Statistics 141 Case Project 

A common problem exists in assessing eye contours as computers have proven to be poor at distinguishing where the eyelid starts and ends. A solution is to try crowdsourcing by asking people to mark an eyelid with a specified amount of points that outline the eyelid’s shape. 

This paper will present an analysis of the accuracy of individuals hired in identifying the shape and location of an eyelid in comparison to an expert’s review. Specifically, analysis will be conducted to answer questions such as accuracy of participants’ responses to the expert’s golden standard, the optimal number of clicks participants should make to achieve accurate shapes, and comparisons between experts’ gold standards of different eyes as well as between the same eye pre versus post operation. Analysis is conducted by comparing polynomials across responses which are created by fitting each individual’s clicks to a 2nd degree polynomial. 

The two-sample Kolmogorov-Smirnov test is used to conclude whether or not any two curves can be treated as the same by testing the null hypothesis of whether the samples (points that construct the fitted polynomial) come from a population from the same distribution. 

The KS test demonstrates that while in most cases, each participant’s clicks are not consistent with the expert’s based on a significance level of 0.05, upon averaging all clicks made by participants and removing any outliers, the polynomial fitted to this averaged matrix of points and the expert’s polynomial come from the same distribution (p-value of 1). As most participant’s individual polynomials are statistically significant and reject the null hypothesis, no conclusion can be made on the optimal number of clicks one needs to make to outline the shape of an eyelid. However, the KS test shows that using 13 or more clicks to outline an eyelid will offer results which are more consistent with the expert’s. Furthermore, the KS test demonstrates that there is an overall difference in shapes of different individuals’ eyes and eyes of the same individual pre and post-operation. 

The next goal of this study is to set precise parameters for participants to follow in order to obtain the most accurate measurements.


Questions Addressed:
1) How well do participants’ responses map to the expert’s golden standard?
2) What is the optimal number of clicks for participants?
3) Are experts’ gold standards always uniform?
4) Is there a significant change in shape in pre- versus post-operation?
