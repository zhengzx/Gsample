# Gsample
Practical Error of convolution Gaussian Sampling

The code is based on NTL library.

It is used to test practical errors of convolution Gaussian sampling with different parameters under different precisions.

Which means x1,x2 are sampled by discrete Gaussian sampling centered around 0 and standard deviations are s1,s2.

And x is sampled by discrete Gaussian sampling centered around 0 and standard deviation is s=sqrt((as1)^2+(bs2)^2).

Let x'=ax1+bx2, we check the differences between the distribution of x' and x and denote it as the practical error of convolution Gaussian sampling. 

More specifically, 

In example 1: a=11,b=1,s1=s2=19.53sqrt(2pi),t=3-8 ,pre=53-200 (experiment about \Delta_{SD} and \Delta_{KL})

In example 2: a=4,b=3,s1=s2=34,t=3-8,pre=53-200 (experiment about \Delta_{RE} and \Delta_{ML})

In example 3: a=11,b=1,s1=s2=19.53sqrt(2pi),t=5.35 ,pre=53-200 ([PDG14]: pre=72, Modified [PDG14]: pre=130)

In example 4: a=4,b=3,s1=s2=34,t=6,pre=53-200 ([MW17]: pre=60, \Delta_{ML}\le 2^{-55}  Modified [MW17]: pre=113, \Delta_{KL}\le 2^{-110})


Other parameters except a,b,s1,s2 are:  t defines the scope for the variable x as [-ts,ts], and pre (precision) means the precision used in calculating.

All examples are independent and can be chosen to test or not by switches, for example:

"EX2 = 1;" means example 2 will be run and "EX3 = 0;" means example 3 will be skipped.

Similarly, other experiments with arbitrary parameters can also be implied.
