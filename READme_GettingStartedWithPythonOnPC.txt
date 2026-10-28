Step-by-step directions to prepare conda env, by LRW on Kyle machine and my laptop
1.	Open Anaconda Prompt from Windows search bar
2.	Paste ONE AT A TIME: 
conda create -n modis312 -c conda-forge python=3.12
conda activate modis312
conda install -c conda-forge earthengine-api matplotlib geopandas shapely pandas
python -m pip install geemap*++++++ and any other packages ++++++*
3.	Follow prompts 
4.	Choose modis:  $ conda activate modis312
5.	Follow prompts 
6.	Close and open Anaconda prompt
7.	Type conda activate modis AND THEN python -c "import ee, geemap, matplotlib, geopandas, shapely, pandas; print('Imports successful')"
8.	Open VS code editor
9.	In the notebook, select the modis312 kernel from the kernel picker near the top right. Choose Select another kernel > Python environments.
10.	Run the import cell again. If the notebook was already open during installation, restart its kernel first.
11.	The notebook also uses Google Earth Engine. After selecting the correct kernel, you may need to run ee.Authenticate() once and complete the browser sign-in. 
a.	Uncomment ee.Authenticate() and log in as uwyo account
b.	Follow prompts
c.	Copy the authentication code (ex 4/1AXlqoi7ZYGrjOXv6Ocj2TKiWxPVD9HmQA9H5MFZacio-OVaCdEh_y0j6750) 
d.	Paste it in the VS code at the top dialogue box
12.	Then run ee.Initialize(project="westwidemuledeer"). Earth Engine access must be enabled for that Google account and project.


Step-by-step directions to produce NDVI data for WWMD for me, by LRW on Kyle machine and my laptop
1.	Open \DroughtSurfing\1_RemoteSensingProcessing_PythonGEE.ipynb into visual studio. May need to adjust the gee authentication info if a different user
a.	If you haven’t already, clone from tlasharr github 
b.	Be sure to start your own branch to keep things from messing with main
c.	For me (LRW), keep using your branch and local file in Documents/Github
2.	In the file adjust the start and end dates to be the date range you want to export to. Right now, we are doing the entire Western US over 9 tiles
3.	These will write to the specified location in My Drive, but that chokes most storage capabilities. So, you will want to pull them onto a local C drive. You can cut and paste the set of 9 for each year at a time, or reactivate the task scheduler task written last time (“Move NDVI tile files every 4 hours” on Kyle PC). This copies from G/NDVI_export to C/Desktop/tmp ndvi/ or similar
4.	Then, begin the R code 1_MergeTilesExportedFromGEE.R which runs in parallel. Ensure that it has the right paths. This will write four of the tifs at a time to a local drive. You will then be able to copy them piecemeal (no more than 2 in a cut, so that it runs parallel). The R code is also running parallel, over 4 cores. This works and doesn’t crush CPU or RAM
5.	To Be continued…

Updating Python
1.	If you created the `modis` environment from the earlier steps, update Python **inside that environment** using **Anaconda Prompt**:
2.	ONE AT A TIME 
conda activate modis
python --version
conda update python
python --version

3.	Review Conda’s proposed changes and enter `y` to confirm. `conda update python` installs the newest Python version that works with the packages already in that environment.

4.	If you want to choose a specific version instead, for example Python 3.12, use:

5.	conda install python=3.12

6.	Afterward, restart the notebook’s kernel in VS Code and make sure it’s still using the `modis` environment. To see your environments, run `conda env list` in Anaconda Prompt.
