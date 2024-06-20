# Fat-Correction
Traditionally, CEST contrast is calculated by the asymmetry analysis, but the presence of fat signals lead to wrong contrast quantification, hence to inaccurate pH measurements. We investigated four post-processing approaches to overcome fat signal influences that enable correct CEST contrast calculation and tumor pH measurements:
- positive method (#1) calculates the contrast considering only the positive part of the Z-spectrum
- linear method (#2) consists in removing the fat frequencies and replacing the missing range with a linear interpolation
- lorentzian method (#3) consists in replacing the negative part of the Z-spectrum by the water pool contribution upon Lorentzian fitting of the spectrum
- interpolation method (#4) corrects the calculated ratiometric values accordingly to the measured fat fraction levels, by interpolating the ratiometric values with cubic splines to correct for the proper pH values

# Getting Started
The project was written in Matlab R2023b (The Math-Works, Inc., Natick, MA, USA). Lorentzian CEST curve fitting was implemented with the open-source Matlab-based code (https://github.com/cest-sources/CEST_EVAL). 


