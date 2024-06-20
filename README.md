# Fat-Correction
Traditionally, CEST contrast is calculated by the asymmetry analysis, but the presence of fat signals lead to wrong contrast quantification, hence to inaccurate pH measurements. We investigated four post-processing approaches to overcome fat signal influences that enable correct CEST contrast calculation and tumor pH measurements:
- positive method (#1) calculates the contrast considering only the positive part of the Z-spectrum
- linear method (#2) consists in removing the fat frequencies and replacing the missing range with a linear interpolation
- lorentzian method (#3) consists in replacing the negative part of the Z-spectrum by the water pool contribution upon Lorentzian fitting of the spectrum
- interpolation method (#4) corrects the calculated ratiometric values accordingly to the measured fat fraction levels, by interpolating the ratiometric values with cubic splines to correct for the proper pH values

# Getting Started
The project was written in Matlab R2023b (The Math-Works, Inc., Natick, MA, USA). Lorentzian CEST curve fitting was implemented with the open-source Matlab-based code (https://github.com/cest-sources/CEST_EVAL). 

## In vitro analysis
Two phantoms were prepared for in vitro validation. The phantoms consisted of a solution containing 30 mM iopamidol, one with the pH adjusted to 6.4 and the other with the pH adjusted to 6.9. Raw Bruker data are provided:
- 20200721_164409_dv_phantom_fat_iopa_3_1_3_pH_6_4.zip:
   * pH = 6.4
   * anatomical folder = 8
   * CEST pre folder = 10
     
- 20210324_094540_dv_iopa_fat_invitro_1_1_1_pH_6_9.zip: 
   * pH = 6.9
   * anatomical folder = 9
   * CEST pre folder = 10
     
In vitro analysis can be carried out by running the script "FatCorrectionVitro.m". Inside the phantom folder, a results folder called "FF" is created which contains the Z-spectra and parametric maps resulting from the analysis of the asymmetric method and the 4 proposed methods.



