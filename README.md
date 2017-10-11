# Gsample
Practical Error of convolution Gaussian Sampling

The code is based on NTL library.
It is used to test practical errors of convolution Gaussian sampling with different parameters under different precisions.
Which means x1,x2 are sampled by discrete Gaussian sampling centered around 0 and standard deviations are s1,s2.
And x is sampled by discrete Gaussian sampling centered around 0 and standard deviation is s=sqrt((ax1)^2+(bx2)^2).
Let x'=ax1+bx2, we check the differences between the distribution of x' and x and denote it as the practical error of convolution Gaussian sampling. 

More specifically, 
In example 1: s1=s2=34, a=b=1
In example 2: s1=s2=19.53, a=11, b=1
In example 3: s1=s2=19.53, a=1,...,11, b=1
In example 4: s1=s2=215.73/7, a=7, b=1
In example 5: s1=s2=34, a=4, b=3
Other parameters are: sp stands for the smoothing parameter, t defines the scope for the variable x as [-ts,ts], and precision means the precision used in calculating.

All examples are independent and can be chosen to test or not by switches, for example:

"EX2 = 1;" means example 2 will be run and "EX3 = 0;" means example 3 will be skipped.

Similarly, other experiments with arbitrary parameters can also be implied.
