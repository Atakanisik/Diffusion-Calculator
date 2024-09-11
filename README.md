# Diffusion-Calculator
Diffusion Calculator is for Calculation IVIM and ADC parameters from DWI-MRI images espacially for SIEMENS vendor DICOM tags.

Loading DICOM Series:
-Use Load DICOM button on bottom right in Calculator TAB for importing DICOM series. You can import multiple DICOM series at once. When Loading is finished you can see Ågdicom installedÅh warning on bottom left.




You can see your DICOM Collection on DICOM TAB

-Calculator:
--You can load DWI trace images on Calculator. If you try load another sequence you will get error but you can load any sequence on the Viewer Directly.

DWI trace sequence loaded to calculator when you press load to calculator button, after that you can see your images on calculator tab with 2 knob (b-Value and Slice Number).


Draw ROI and Multiple Draw buttons for image region selection and you can see ROI types right corner on  dropdown menu. Draw ROI is for single ROI selection. For Multiple ROI you can use draw Multi ROI button. After selecting multiple ROI to interrupt drawing you can switch the button below of the Multi ROI button ÅgoffÅh and then click the most bottom right corner image to stop drawing ROIs.


A window will be appeared when you complete the selection and then immediately you can see signal decay graphics on the right panel.



Then You can calculate IVIM and ADC parameters via just pressing the calculate button.
2 Algorithm is included for IVIM free-fit and Segmented-fit. You can select the algorithms from radio buttons. After Calculation will be completed you can see all parameters and fitting curve.If you use all algorithms you can generate maps for specific regions.



Pressing Generate Maps button will generate IVIM and ADC maps for all algorithms. Be sure to use single or multi ROI


After calculation is done you can see your maps in Maps TAB.


If you have support packages you can use segmenter with SAM model.
Also Application include BOLD(R2*) Processor for recalculate R2* images from multi TE Series you can monitor the oxygenization from these Series.

You can use Viewer for all sequences.

 if you use this prototype application please refer me as:
Research Assistant Atakan I??k
Baskent University Biomedical Engineering Department.
--ORC-ID:0000-0001-5433-4442