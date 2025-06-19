A MATLAB-based code for Zld-sfGFP and GAF-sfGFP image analysis and fitting models to data for the paper:

Dima, S. S. & Reeves, G. T. (2024). Bulk-level maps of pioneer factor binding dynamics during Drosophila maternal-to-zygotic transition Development, dev.204460. DOI: 10.1242/dev.204460

Installation instructions for the RICS analysis:

- This code requires the BioFormats plugin for Matlab, which can be accessed from https://www.openmicroscopy.org/bio-formats/. 
- This pipeline was created using BioFormats version 6.10.1. Updated versions of BioFormats may change functionality, and have not been tested.
- Install BioFormats on your computer first. 
- Download all contents into a single folder on your local computer. 
- Open MATLAB. 
- Click "Home>Set Path" on the top bar of the window. This will cause a "Set Path" window to appear. 
- Click the "Add With Subfolders" option.
- You will need to download all the code in this folder of this repository and make it available to Matlab on your computer
- Note that a lot of the code in this folder can also be found at: https://github.com/gtreeves/RICS_timecourse_pipeline
- Additionally, this pipeline uses functionality found in the RICS pipeline (non-timecourse), which can be found at: https://github.com/gtreeves/RICS_pipeline
- If you wish to run the RICS_timecourse pipeline on the images from the paper, you will need to download the image data set from the Texas Data Repository
- The folder you place the image data sets in should be reflected in the "Tiles.xlsx" excel file, so you will have to edit Column A to achieve that
- Next, run "script_analyze_Zld" or "script_analyze_GAF"
- The results should be saved in a mat file. If done correctly, the data in these mat files should reflect the mat files already stored in the "Generate figures" folder





