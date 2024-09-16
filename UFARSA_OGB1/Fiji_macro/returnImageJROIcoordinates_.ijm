// This macro returns the x-y coordinates of all pixels for each ROI in MATLAB notation.
// All types of ImageJ ROIs are supported.
// Tested with Fiji 1.53c and Java version "1.8.0_333", Java(TM) SE Runtime Environment (build 1.8.0_333-b02),
// Java HotSpot(TM) 64-Bit Server VM (build 25.333-b02, mixed mode)
// Author: Dr. Knut Kirmse, 2018-2022.
// ------------------------------------------------------------
// Clean-up
if (isOpen("ROI Manager"))
	{
    selectWindow("ROI Manager");
    run("Close");
    }
run("Close All");

// --- File selection: ROI-manager ROI list ('*.roi' oder '*.zip') ---
width = 1024; // Note: ImageJ requires an appropriately sized image to extract coordinates from ROIs
height = 1024;
Dialog.create("Image dimensions")
Dialog.addMessage("The given dimensions must be equal to or larger than\nthose of the image in which the ROIs were selected.");
Dialog.addNumber("Width [pixels]: ", width)
Dialog.addNumber("Height [pixels]: ", height)
Dialog.show();
width = Dialog.getNumber();
height = Dialog.getNumber();
ROIpath = File.openDialog("Select a ROI list file");
open(ROIpath);

// if there is only one ROI, this needs to be added to ROI Manager first
nROIs = roiManager("count");
if (nROIs == 0)
	{
	roiManager("Add");
	nROIs = roiManager("count");
	}
run("Close All");
roiManager("Associate", "false");
roiManager("Centered", "false");
roiManager("UseNames", "false");
newImage("Untitled", "8-bit white", width, height, 1);
rename("Placeholder");

// determine maximum pixel count per ROI
maxNrows = 0;
for (k = 0; k < nROIs; k++)
	{
	roiManager("Select", k);
	Roi.getContainedPoints(xpoints, ypoints);
	if (maxNrows < xpoints.length)
		{
		maxNrows = xpoints.length;
		}
	}

// fill table with coordinates and fill empty cells with NaNs
Table.create("X-Y coordinates of ROIs");
for (k = 0; k < nROIs; k++)
	{
	roiManager("Select", k);
	Roi.getContainedPoints(xpoints, ypoints);
	for (m = 0; m < xpoints.length; m++)
	{
		xpoints[m] = xpoints[m] + 1; // in Matlab, the first element has index 1
		ypoints[m] = ypoints[m] + 1; // in Matlab, the first element has index 1
	}
	aux = newArray(maxNrows - xpoints.length);
	Array.fill(aux, NaN);
	rowvalues = Array.concat(ypoints, aux);
	colvalues = Array.concat(xpoints, aux);
	Table.setColumn("row_ROI " + k+1, rowvalues);
	Table.setColumn("col_ROI " + k+1, colvalues);
	}

// Save x-y coordinates as a TXT file in the same directory
Table.showRowNumbers(false);
Table.showRowIndexes(false);
Table.saveColumnHeader(true);
Table.save(ROIpath + ".txt");

// Close
if (isOpen("Placeholder"))
	{
    selectWindow("Placeholder");
    run("Close");
    }