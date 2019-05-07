Matlab code for Gaussian Processes Probabilistic Amplitude
Demodulation (GPPAD).

Copyright Richard E. Turner, University of Cambridge, 2010

This approach to demodulation is described in: 

Turner and Sahani 2010, 'Demodulation as an Inference Problem', TASLP.
Turner 2010, 'Statistical models for natural sounds', UCL PhD thesis. 
Turner and Sahani 'Probabilistic Amplitude Demodulation' 2007. 

These papers can be found in the publications directory.

The main function is GPPAD.m so, for additional information, see the
help (i.e. type 'help GPPAD' into the matlab command line).

In addition to this help, demonstrations of some of the options can be
found in the \Demos directory:

GPPADDemo1.m - demodulates a one-dimensional signal using GPPAD. Some
parameters are learned from the signal, but not the time-scale of the
modulation which must be specified by the user. This is a good place
to get started - a (relatively) fast and simple setup.

GPPADDemo2.m - demodulates a one dimensional signal using GPPAD and
learns the time-scale of the modulation from the signal together with
the other parameters. This example also shows how to request
error-bars in the modulators to be returned.
 
GPPADDemo3.m - uses parameters that have been derived previously on a
clean signal, to denoise and fill in parts of a new noise corrupted
signal. Error-bars are returned again for both the original and noisy
estimated modulator.

GPPADDemo4.m - demodulates a one dimensional signal using GPPAD and
learns the time-scale of the modulation and the marginal parameters
like GPPADDemo2.m. However, the time-scales are learned using a slow
grid based method here. Although the simulation is computationally
more intensive, the advantage is that it also returns the likelihood
of the time-scales over a range of user-specified values. In the
speech example used in the demo, the likelihood is multimodal, which
indicates that there are two possible ways to demodulate the signal;
one at the rate of the pitch, and one at the rate of the phonemes.

GPPADDemo5.m - demonstrates sub-band demodulation. A sentence sound is
passed through a gammatone filter bank, and then GPPAD is used to
demodulate each filter output. The user specifies the time-scale in
this example, but it can also be learned (see GPPADDemo2.m). The other
parameters are learned for each filter.

GPPADDemo6.m - demonstrates the cascade model, that is recursive
demodulation at longer and longer time-scales. Here a spoken sound is
demoduled at the rate of the pitch, the phonemes, and the sentences.
