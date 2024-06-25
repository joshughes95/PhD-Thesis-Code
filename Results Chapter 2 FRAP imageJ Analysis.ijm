dir1 = getDirectory("Choose source directory ");
list = getFileList(dir1);
dir2 = getDirectory("Choose destination directory ");

Experiment = 3

imageId = getImageID();
//Loop_over_all_ROIs
number_of_rois = roiManager("count");
run("Set Measurements...", "mean redirect=None decimal=3");
setOption("ShowRowNumbers",false);
setOption("ShowRowIndexes", false);
for (i=0; i<number_of_rois+1; i++) {
//Runs_ZProfile_for_all_ROIs_adds_to_results
	selectImage(imageId); //need to reselect image in each cycle
	if (i<number_of_rois) {
		roiManager("Select", i);
		run("Plot Z-axis Profile");
	} else {
		run("Select None");
		selectImage(imageId);
		run("Duplicate...", "duplicate");
		setOption("BlackBackground", true);
		run("Make Binary", "method=Huang background=Dark calculate black");
		run("Analyze Particles...", "size=20-Infinity pixel add slice");
		close();
		roiManager("select", Array.slice(Array.getSequence(roiManager("count")),number_of_rois,roiManager("count")));
		roiManager("combine");
		roiManager("add");
		roiManager("select", Array.slice(Array.getSequence(roiManager("count")),number_of_rois,roiManager("count")-1));
		roiManager("delete")
		selectImage(imageId);
		roiManager("select", roiManager("count")-1);
		run("Plot Z-axis Profile");
	}
	xpoints = newArray ();
	ypoints = newArray ();
	Plot.getValues(xpoints, ypoints);
//Move data from plot values into results
	ROI_name = "ROI_"+i;
	for (j = 0; j < xpoints.length; j++) {
		setResult ("Time", j, xpoints[j]);
		setResult (ROI_name, j, ypoints[j]);
	}
	updateResults();
	close();
	
saveAs ("Results", dir2+"FRAP_CSV_"+Experiment+".csv");
}

//Measure Area of ROIs and whole Cells and Save to seperate CSV
run("Clear Results");
run("Set Measurements...", "area redirect=None decimal=3");
selectImage(imageId);
roiManager("Select", 0);
run("Measure");
roiManager("Select", roiManager("count")-1);
run("Measure");
saveAs ("Results", dir2+"AREA_CSV_"+Experiment+".csv");
run("Clear Results");