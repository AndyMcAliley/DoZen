# README #

### What is this repository for? ###

DoZen (pronounced "do zen") is for processing, visualizing, and exploring electromagnetic data that is stored in Zonge's .z3d format. It was created specifically for processing time-lapse controlled source electromagnetic data.

### How do I get set up? ###

## Get the code:
```
git clone https://github.com/AndyMcAliley/DoZen.git
```

## Install:
First install DoZen and all necessary dependencies. I recommend using Anaconda. If you do, you can create a new environment just for DoZen:
```
conda env create -f environment.yml
conda activate dozen
```
If you want interactive plots in Jupyterlab, you'll need extensions.
See https://jupyterlab.readthedocs.io/en/stable/user/extensions.html for details.
Run the following commands to install the necessary conda packages for extensions:
```
conda install -c conda-forge ipympl
conda install -c conda-forge nodejs
```
Next, run
```
jupyter lab
```
From Jupyter lab, you can enable extensions by clicking on the commands pane on the left, then searching for Extension Manager.
Once extensions are enabled, run the following commands to install the necessary Jupyterlab extensions:
```
jupyter labextension install @jupyter-widgets/jupyterlab-manager
jupyter labextension install jupyter-matplotlib
jupyter labextension install jupyterlab_bokeh
jupyter labextension install @pyviz/jupyterlab_pyviz
```

### How do I use it? ###

## From within Python:
First and foremost, you can import and use the library from within Python:
```
import dozen
```

To read a z3d file:
```
dozen.z3dio.read_z3d(filename)
```

To read metadata for all files in a directory into a dataframe:
```
dozen.z3d_directory.directory_info()
```
OR
```
dozen.z3d_directory.directory_info(initialdir=path_to_directory,ask_dir=False)
```

## From a script:
In addition, you can run one of the scripts which open dashboards:

To view a single z3d timeseries:
Open a terminal in the DoZen folder, and then
```
cd ipynb
panel serve --show Plot_z3d.ipynb
```
Click the "Load z3d" button

To view a map of all file locations,
Open a terminal in the DoZen folder, and then
```
cd ipynb
jupyter lab
```
Next, change the file locations and other campaign-specific variables.
Exit Jupyter, and then
```
panel serve --show CSEM_Preprocess.ipynb
```
From here you can select a group of locations and save its average.
Once you've averaged several stations, you can export those locations as a csv.
You can also export all the receiver z3d file metadata in one large csv.

To view an interactive timeline of CSEM data in a browser window:
edit scripts/timeline_settings.csv to point to the relevant directories and files.
Open a terminal in the DoZen folder, then
```
cd scripts
panel serve --show timeline.py
```
Select a campaign from the dropdown menu.
Click on the timeline to select a transmission.
Click the plot button to plot the selected transmission.

To process data:
Write a script (you can pattern it after process_example.py)
Then,
```
python your_processing_script.py
```

To rotate data:
Open a terminal in the DoZen folder, then
```
cd scripts
python rotate.py
```


You can also run these as jupyter notebooks (although some of the interactivity like time series plotting when you click a part of a timeline, may not work in a notebook).
