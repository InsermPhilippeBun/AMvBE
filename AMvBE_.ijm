// Macro version written by Philippe BUN (IPNP, INSERM U1266)
// Minimal requirements - ImageJ 1.51u
// TurboReg pluging from EPFL and MultiStackReg from Brad Brusse

macro "AMvBE" {

// ---- No saved ROIs prior to recording positions 
if (isOpen("ROI Manager")) { 
	roiManager("reset");
	roiManager("Show All with labels");
} else {
	roiManager("show none");
	roiManager("Show All with labels");
}

// ---- For the table result
run("Clear Results");
// ---- For the log window
print("\\Clear");
// ---- For the images
run("Close All");


// ----- Directory -----

open();
name_file = File.nameWithoutExtension;
full_name_file = File.name;
dir = getDirectory("Saving directory");
save_folder = dir + name_file + "_results/";
if ( File.exists(save_folder) < 1 ) File.makeDirectory( save_folder );

// ----- Settings -----

setForegroundColor(0, 0, 0);
setBackgroundColor(255, 255, 255);
print("\\Clear");
run("Clear Results");
run("Set Measurements...", "area mean min perimeter integrated display redirect=None decimal=5");
run("Plots...", "width=530 height=300 font=12 draw draw_ticks auto-close minimum=0 maximum=0 interpolate sub-pixel");
run("Set... ", "zoom=150"); 

// ----- User input box -----
Dialog.create("Exocytosis event -manual selection -automatic confirmation");
Dialog.addCheckbox("Bleaching correction ?", false);
Dialog.addCheckbox("Registration to correct cell movement ?", true);
Dialog.addCheckbox("Have you some pre-saved ROIs ?", false);
Dialog.addCheckbox("Generation of exocytosis event stacks", true);
Dialog.addNumber("Pixel size in micrometer", 0.16); //pixel size
Dialog.addNumber("Frame rate in second", 0.50); //frame rate
Dialog.addNumber("Time window for event correction (in nb of frames)", 6); // Number of frames to spot the exocytosis event
Dialog.addNumber("Condition for fluorescence decrease after exocytosis (min. nb of frames)",3);
Dialog.addNumber("Condition for fluorescence increase prior to exocytosis (max. nb of frames)",3);
Dialog.addNumber("Min. event size to be considered (in micrometer)",0.4);
Dialog.addNumber("If the decay condition is not fulfilled, \r threshold value (N) for exocytosis event rescue such that Mean + N*STD", 3);
Dialog.addNumber("Manual threshold for vesicle movement fixed by users (unit in number of vesicle diameter)", 2);

Dialog.show();
Bleaching = Dialog.getCheckbox();
Registration = Dialog.getCheckbox();
saveROI = Dialog.getCheckbox();
Display = Dialog.getCheckbox();
pixel_size = Dialog.getNumber();
frame_rate = Dialog.getNumber();
duration_test = Dialog.getNumber();
decrease_condition = Dialog.getNumber();
increase_condition = Dialog.getNumber();
size_condition = Dialog.getNumber();
I_threshold = Dialog.getNumber();
manual_threshold = Dialog.getNumber();

// ----- END -----

//----- Settings -----
ImageID = getImageID(); 
run("Select None");
run("Line Width...", "line=1");
run("Options...", "iterations=1 count=1 edm=Overwrite");
number_slices = nSlices();

// ----- Normalization -----
//r_correction = 2*size_condition;
r_correction = size_condition;

//half_length = 10; // in pixel for plot profile == FWHM
half_length = 4*size_condition/pixel_size;
angular_resolution = 30; //in degrees

// ----- Arrays -----
Rsquared_array = newArray(360/angular_resolution);
counter_array = newArray(360/angular_resolution);
Rsquared_array_sampled = newArray(360/angular_resolution);
counter_array_sampled = newArray(360/angular_resolution);

// ----- Bleaching correction -----
if (Bleaching == true) {
waitForUser("Bleach correction","Post-processing - Bleach correction.");
Image_preprocessing = getImageID();
run("Bleach Correction", "correction=[Exponential Fit]");
ImageID = getImageID();
selectImage(Image_preprocessing); close();
print("\\Clear");
}

// ----- Registration -----
if (Registration == true) {
	selectImage(ImageID);
	List.setCommands;
    if (List.get("TurboReg ")!="") {
    	if(List.get("MultiStackReg")!="" ) {
           	setSlice(1); run("MultiStackReg", "stack_1=full_name_file action_1=Align file_1=[] stack_2=None action_2=Ignore file_2=[] transformation=Translation ");
			print("\\Clear");
    	} else {
    		exit("MultiStackReg plugin was not (properly) installed. Run will be terminated here.");
    	}
    } else {
    	exit("TurboReg plugin was not (properly) installed. Run will be terminated here.");
    }
}

// ---- summary ---
print("");
print("User-defined parameters");
print("-----------------------");
if (Bleaching) { print("Bleaching correction: YES"); } else { print("Bleaching correction: NO"); }
if (Registration) { print("Registration: YES"); } else { print("Registration: NO"); };
print("Pixel size is: " +pixel_size+ " micrometer");
print("Frame rate is: " +frame_rate+ " s");
print("Size of the window for correction in time to spot the exocytosis event is: " +duration_test*frame_rate+ " s");
print("Size of the window for correction in space to spot the exocytosis event is: " +r_correction+ " micrometer");
print("You choose an angular resolution of: "+angular_resolution+ " to precisely determine the size of vesicles");
print("");
print("User-defined conditions");
print("-----------------------");
print("The size threshold is: "+size_condition+ " micrometer");
print("The vesicle should be confined within a circular surface of " +manual_threshold+ " times the vesicle diameter in micrometers"); 
print("The fluorescence increase/burst should last less than "+increase_condition+ " frames"); 
print("The fluorescence decrease should last at least "+decrease_condition+ " frames"); 
print("The value for fluorescence rescue if the fluorescence decrease condition was not fulfilled was fixed to the average intensity + "+I_threshold+ " * the standard deviation");
print("----------------");

//--------------------------------------------------------------------------------

if (saveROI == false) {
waitForUser("Manual selection","Pinpoint the exocytosis event (Left click).\nRight click to exit.");

selectImage(ImageID);
// ----- Event selection -----

//setOption("DisablePopupMenu", true);
setTool("rectangle"); setTool("point");
leftButton=16; 
rightButton=4;
x2=-1; y2=-1; z2=-1; flags2=-1; 
getCursorLoc(x, y, z, flags); 
Outclick = false; 
while (flags&rightButton==0){ 
	getCursorLoc(x, y, z, flags);      
	if (flags&leftButton!=0) { // Wait for it to be released 
		Outclick = true; 
        } else if (Outclick) {
        Outclick = false;
        if (x!=x2 || y!=y2 || z!=z2 || flags!=flags2) {
        	makePoint(x,y);
            roiManager("Add"); 
            setOption("DisablePopupMenu", true);
        }
	}
}                           
roiManager("Show All without labels");
roiManager("Show None");
//--------------------------------------------------------------------------------

run("Select None");
roiManager("save", save_folder + name_file +"_before correction" + ".zip");



} else {
 	waitForUser("Pre-chosen ROIs","Point towards your ROIs after pressing OK");
 	open();
 	run("Select None");
}

roiManager("Show All without labels");
roiManager("Show None");
waitForUser("Event analysis","Each event will be analyzed.");

// ----- Settings -----
setBatchMode(true);
number_ROI = roiManager("Count"); roiManager("Deselect");
selectImage(ImageID); 
run("Select None");
run("Set Scale...", "distance=1 known=1 pixel=1 unit=unit global");
run("Set Scale...", "distance=0");
Stack.getDimensions(width, height, channels, slices, frames);
run("Properties...", "channels="+channels+ " slices="+slices+ " frames="+frames+ " unit=unit pixel_width="+pixel_size+ " pixel_height="+pixel_size+ " voxel_depth="+frame_rate+ " frame="+frame_rate);


for (i=0; i<number_ROI; i++) {
	setResult("ROI", i, i+1);
	selectImage(ImageID);
	roiManager("Select",i);
	user_timeposition = getSliceNumber(); 
	getSelectionCoordinates(x_user,y_user);
	setResult("user X",i,x_user[0]);
	setResult("user Y",i,y_user[0]);
	
	// Centering correction -max intensity
	makeOval(x_user[0]-r_correction/pixel_size, y_user[0]-r_correction/pixel_size, 2*r_correction/pixel_size, 2*r_correction/pixel_size);
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	run("Find Maxima...", "noise="+std+" output=[Point Selection]");
	getSelectionCoordinates(corr_x, corr_y);
	delta_r = sqrt(pow(corr_y[0]-y_user[0],2)+pow(corr_x[0]-x_user[0],2))*pixel_size;
	
	setResult("corr X",i,corr_x[0]);
	setResult("corr Y",i,corr_y[0]);
	setResult("delta R",i,delta_r);

	//defining vesicle size (FWHM)
	loop = -1;
	counter = 0;
	makeLine(corr_x[0]-half_length,corr_y[0],corr_x[0]+half_length,corr_y[0],1);
	
	setBatchMode(true);
	do {
		selectImage(ImageID);
		loop = loop + 1;
		counter = counter + angular_resolution;
		run("Rotate...", "  angle=counter");
		run("Plot Profile");
		profile_rotate = getImageID();
		Plot.getValues(x, y); close();
		R2 = FitOutcome(x, y, "Gaussian", "R2");
		Rsquared_array [loop] = R2;
		counter_array [loop] = counter;
	} while (counter < 360);
	Array.getStatistics(Rsquared_array,min,max);

if (max > 0.8){ //condition to find the best fit for vesicle size -here fixed to 0.9

	index_angle = IndexMax(Rsquared_array,max);
	optimal_angle = angular_resolution*(index_angle[index_angle.length-1]+1);
			
	selectImage(ImageID);
	run("Select None");
	makeLine(corr_x[0]-half_length,corr_y[0],corr_x[0]+half_length,corr_y[0],1);
	
	run("Rotate...", "  angle=optimal_angle");
	run("Plot Profile");
	Plot.getValues(x, y); close();
	FWHM = FitOutcome(x, y, "Gaussian", "FWHM");

	//Defining vesicle size with proper sampling for fits
	selectImage(ImageID);
	run("Select None");
	sampling_length = 2*round(FWHM/pixel_size);
	makeLine(corr_x[0]-sampling_length,corr_y[0],corr_x[0]+sampling_length,corr_y[0],1); //adjusted length
	run("Rotate...", "  angle=optimal_angle"); // adjusted orientation
	run("Plot Profile");
	Plot.getValues(x, y);
	close();
	FWHM = FitOutcome(x, y, "Gaussian", "FWHM");
	R2 = FitOutcome(x, y, "Gaussian", "R2");
	adjR2 = 1 - ((1-R2)*(x.length-1)/(x.length-3-1));

	//optimization
	loop = -1;
	counter = 0;
	makeLine(corr_x[0]-sampling_length,corr_y[0],corr_x[0]+sampling_length,corr_y[0],1);
	setBatchMode(true);
	do {
		selectImage(ImageID);
		loop = loop + 1;
		counter = counter + angular_resolution;
		run("Rotate...", "  angle=counter");
		run("Plot Profile");
		profile_rotate = getImageID();
		Plot.getValues(x, y);
		close();
		R2 = FitOutcome(x, y, "Gaussian", "R2");
		Rsquared_array_sampled [loop] = R2;
		counter_array_sampled [loop] = counter;
		} while (counter < 360);
	setBatchMode(false);	
	Array.getStatistics(Rsquared_array_sampled,min,max);
	index_angle_sampled = IndexMax(Rsquared_array_sampled,max);
	optimal_angle_sampled = angular_resolution*(index_angle_sampled[index_angle_sampled.length-1]+1);

	selectImage(ImageID);
	run("Select None");
	makeLine(corr_x[0]-sampling_length,corr_y[0],corr_x[0]+sampling_length,corr_y[0],1);
	run("Rotate...", "  angle=optimal_angle_sampled");
	run("Plot Profile");
	Plot.getValues(x, y); close();
	FWHM = FitOutcome(x, y, "Gaussian", "FWHM");

	// Correction in time
	selectImage(ImageID);
	makeOval(corr_x[0]-FWHM/pixel_size, corr_y[0]-FWHM/pixel_size, 2*FWHM/pixel_size, 2*FWHM/pixel_size);
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	tolerance = max-min;
	start_z = user_timeposition-duration_test;
	end_z = user_timeposition+duration_test;
	 
	if (start_z <= 0) {
		start_z = 1;
	}
	if (end_z > nSlices()) {
		end_z = nSlices();
	}
	
	run("Duplicate...", "duplicate range=start_z-end_z");
	Zprofile_image = getImageID(); 
	run("Plot Z-axis Profile"); 
	Zprofile = getImageID();
	Plot.getValues(x,y);
	close();
	array_zprofile = newArray(x.length);
	for (k=0; k<x.length; k++) {
        array_zprofile[k] = y[k];
    }
	correction = Array.findMaxima(array_zprofile,tolerance/10);
	
if (correction.length >= 1) { // additional condition to prevent empty correction
	
	
	if (start_z <= 0) {
		corr_timeposition = correction[0]+start_z;
	}
	if (start_z > 0) {
		start_z = user_timeposition-duration_test;
		corr_timeposition = correction[0]+start_z;
	}
	selectImage(Zprofile_image); close();
	
	
	setResult("user Slices",i,user_timeposition);
	setResult("corr Slices",i,corr_timeposition);

// Higher condition on the time shift to prevent detection of later events on the same (+/-r_correction) location
time_shift = abs(corr_timeposition-user_timeposition); //in frames
if (time_shift <= 3*decrease_condition) { // the condition relies on the decrease_condition value

	// Re-correction in space
	if (user_timeposition != corr_timeposition) {
		roiManager("Select",i);
		setSlice(corr_timeposition);
		makePoint(corr_x[0],corr_y[0]);
		roiManager("Update");
		makeLine(corr_x[0]-sampling_length,corr_y[0],corr_x[0]+sampling_length,corr_y[0],1);
		run("Set Scale...", "distance=1 known=pixel_size pixel=1 unit=unit global");
		run("Plot Profile");
		Plot.getValues(x, y); close();
		Fit.doFit("Gaussian", x, y); 
		offset = Fit.p(0);
		amp_offset = Fit.p(1);
		mean = Fit.p(2);
		std_deviation = Fit.p(3);
		tolerance = amp_offset-offset;  
		R2 = Fit.rSquared; 
		adjR2 = 1 - ((1-R2)*(x.length-1)/(x.length-3-1));
		FWHM = 2*pow(2*log(2),0.5)*std_deviation; // 
		//close();			
	} else {
		corr_timeposition = user_timeposition;
		roiManager("Select",i);
		setSlice(corr_timeposition);
		makePoint(corr_x[0],corr_y[0]);
		roiManager("Update");
	}
	setResult("Size-FWHM",i,FWHM);

	// Duration of the event -mean lifetime -exponential decay
	selectImage(ImageID);
	setSlice(corr_timeposition);
	makeOval(corr_x[0]-FWHM/pixel_size, corr_y[0]-FWHM/pixel_size, 2*FWHM/pixel_size, 2*FWHM/pixel_size);
	min_dup = corr_timeposition;
	max_dup = min_dup + duration_test; //duration_test as step 0
	run("Duplicate...", "duplicate range=min_dup-max_dup");
	Decay = getImageID();
	selectImage(Decay); rename("Decay");  setSlice(nSlices);
	
if (nSlices > 1) { // Condition on the number of slices to prevent bugs
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	initial_amp_offset = min; 
		
	run("Plot Z-axis Profile");
	Plot.getValues(x, y); close();
	decay = "y = c + a*exp(-(b*x))";
	initialGuesses = newArray(initial_amp_offset, 1, initial_amp_offset);
	Fit.doFit(decay, x, y, initialGuesses);
	decay_rate = Fit.p(1);
		
	mean_lifetime = 1/decay_rate; 

	// Duration of the event -mean lifetime -exponential decay - sampling for better fit

	decay_sampling = 4*mean_lifetime/frame_rate; //in number of frames
	
	if (decay_sampling > 1) {
		max_decay_sampling = round(min_dup + decay_sampling);
	
		selectImage(ImageID);
		run("Duplicate...", "duplicate range=min_dup-max_decay_sampling");
		sampled_dup_decay = getImageID(); rename("Sampled decay");
		selectImage(sampled_dup_decay); setSlice(nSlices);
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		max_fluorescence = max; mean_fluorescence = mean; std_level = std;	
		threshold_fluorescence = mean_fluorescence+I_threshold*std_level;
	
		setResult("Max Fluorescence",i,max_fluorescence);
		setResult("Fluorescence threshold",i,threshold_fluorescence);
	
		// Fitting with proper sampling	
		initial_amp_offset = min; 
		run("Plot Z-axis Profile");
		Zprofile = getImageID();
		Plot.getValues(x, y);close();
		decay = "y = c + a*exp(-(b*x))";
		initialGuesses = newArray(initial_amp_offset, 1, initial_amp_offset);
		Fit.doFit(decay, x, y, initialGuesses);
		sampled_decay_rate = Fit.p(1);
	
		sampled_mean_lifetime = 1/sampled_decay_rate;
		raw_sampled_mean_lifetime = sampled_mean_lifetime/frame_rate;
		selectImage(sampled_dup_decay);

		if (sampled_mean_lifetime/frame_rate > 1 && sampled_mean_lifetime/frame_rate <= number_slices) {
			setResult("decay -raw mean lifetime (frames)", i, raw_sampled_mean_lifetime);
			setResult("decay -mean lifetime (s)",i,sampled_mean_lifetime);
			} else {
			raw_mean_lifetime = mean_lifetime/frame_rate;
			setResult("decay -raw mean lifetime (frames)", i, raw_mean_lifetime);
			setResult("decay -mean lifetime (s)",i, mean_lifetime);
		}
		selectImage(sampled_dup_decay); close();
	} else {
		raw_mean_lifetime = mean_lifetime/frame_rate;
		raw_sampled_mean_lifetime = raw_mean_lifetime;
		setResult("decay -raw mean lifetime (frames)", i, raw_mean_lifetime);
		setResult("decay -mean lifetime (s)",i, mean_lifetime);		
	}
	
	// Fluorescent increase prior to exocytosis event -number of timepoints till max intensity
	selectImage(ImageID);
	run("Reverse");
	min_dup_rev = nSlices() - (corr_timeposition-1); 
	max_dup_rev = min_dup_rev + duration_test; //duration_test as step 0
	run("Duplicate...", "duplicate range=min_dup_rev-max_dup_rev");
	dup_decay_rev = getImageID(); rename("reverse");
	selectImage(dup_decay_rev); setSlice(nSlices);
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	initial_amp_offset_rev = min;
	run("Plot Z-axis Profile");
	decay_plot_rev = getImageID();
	Plot.getValues(x, y); close();
	decay = "y = c + a*exp(b*x)";
	initialGuesses_rev = newArray(initial_amp_offset_rev, 1, initial_amp_offset_rev);
	Fit.doFit(decay, x, y, initialGuesses_rev);
	increase_rate = Fit.p(1);
	mean_lifetime_increase = -1/(increase_rate);
	

	// Duration of to reach the max -mean lifetime -reverse stack exponential decay - sampling for better fit
	
	increase_sampling = 8*(mean_lifetime_increase/frame_rate);
	if (increase_sampling < 4) {
		increase_sampling = 4;
	}
	max_increase_sampling = round(min_dup_rev + increase_sampling);
	
	selectImage(ImageID);
	run("Duplicate...", "duplicate range=min_dup_rev-max_increase_sampling");
	sampled_dup_decay_rev = getImageID(); rename("sampled reverse");
	selectImage(sampled_dup_decay_rev); setSlice(nSlices);
	getRawStatistics(nPixels, mean, min, max, std, histogram);

	
	// Fitting with proper sampling -increase rate	
	initial_amp_offset_rev = min;
	run("Plot Z-axis Profile");
	Plot.getValues(x, y); close();
	decay = "y = c + a*exp(b*x)";
	initialGuesses_rev = newArray(initial_amp_offset_rev, 1, initial_amp_offset_rev);
	Fit.doFit(decay, x, y, initialGuesses_rev);
	sampled_increase_rate = -Fit.p(1);
	
	sampled_mean_lifetime_increase = 1/sampled_increase_rate; 
	raw_sampled_mean_lifetime_increase = sampled_mean_lifetime_increase/frame_rate;
	
	selectImage(sampled_dup_decay_rev); close();
	
	if (sampled_mean_lifetime_increase/frame_rate > 1 && sampled_mean_lifetime_increase/frame_rate <= number_slices) {
		setResult("increase -raw mean lifetime(frames)", i, raw_sampled_mean_lifetime_increase);
		setResult("increase -mean lifetime (s)",i,sampled_mean_lifetime_increase);
		} else {
		raw_mean_lifetime_increase = mean_lifetime_increase/frame_rate;
		setResult("increase -raw mean lifetime(frames)", i, raw_mean_lifetime_increase);
		setResult("increase -mean lifetime (s)",i, mean_lifetime_increase);		
	}

	updateResults();
	
	selectImage(Decay); close();
	selectImage(dup_decay_rev); close();
	

	// Movement of the vesicles right after the burst -max intensity peak in time 
	selectImage(ImageID);
	run("Reverse");
	size_tracking_movement = round(raw_sampled_mean_lifetime); 
	condition_image_size = corr_timeposition + size_tracking_movement;


if (condition_image_size <= nSlices && size_tracking_movement > 0) {	// conditions to avoid bug because of not enough timepoint for movement tracking
	mov_x_array = newArray(size_tracking_movement); mov_x_array [0] = corr_x[0];
	mov_y_array = newArray(size_tracking_movement); mov_y_array [0] = corr_y[0];
	pix_array = newArray(size_tracking_movement);
	
	for (j=1; j < size_tracking_movement; j++) {
		setSlice(corr_timeposition+j);
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		
		run("Find Maxima...", "noise="+std+" output=[Point Selection]");
		getSelectionCoordinates(mov_x, mov_y);
		mov_x_array [j] = mov_x[0];
		mov_y_array [j] = mov_y[0];
		pix_array [j] = getPixel(mov_x[0],mov_y[0]);
		run("Select None");
		makeOval(mov_x[0]-FWHM/pixel_size, mov_y[0]-FWHM/pixel_size, 2*FWHM/pixel_size, 2*FWHM/pixel_size);
	}
	distance = newArray(size_tracking_movement);
	distance [0] = 0;
	for (r = 1; r < size_tracking_movement; r++) {
		distance [r] = sqrt(pow(mov_x_array[r] - mov_x_array[0],2) + pow(mov_y_array[r] - mov_y_array[0],2))*pixel_size;
	}
	max_distance = maxOfArray(distance);
	
	if (max_distance > manual_threshold*FWHM){
			setResult("Flag",i,"moving");
			setResult("Distance", i, max_distance);
			updateResults();
	} else {
			setResult("Flag",i,"still");
			setResult("Distance", i, max_distance);
			updateResults();
	}


	
}  else {
	roiManager("Select",i);
	print("ROI_"+i+1+"\r was not analyzed due to a decay either way too short (< 1 frame) or too long (more than "+number_slices+ " frames) to reliably track the event");

	setResult("Flag",i,"not analyzed");
	updateResults();
}
} else {
	roiManager("Select",i);
	print("ROI_"+i+1+"\r was not analyzed due to a fluorescence intensity decay faster than 1 image!");

	setResult("decay -raw mean lifetime (frames)",i,NaN);
	setResult("decay -mean lifetime (s)",i,NaN);
	setResult("increase -raw mean lifetime(frames)",i,NaN);
	setResult("increase -mean lifetime (s)",i,NaN);
	setResult("Flag",i,"not analyzed");
	updateResults();
}

} else {
	roiManager("Select",i);
	print("ROI_"+i+1+"\r was not analyzed. It could be the same event on the same location pinpointed twice! The time shift is "+time_shift*frame_rate+" s! ");

	setResult("Size-FWHM",i,FWHM);
	setResult("decay -raw mean lifetime (frames)",i,NaN);
	setResult("decay -mean lifetime (s)",i,NaN);
	setResult("increase -raw mean lifetime(frames)",i,NaN);
	setResult("increase -mean lifetime (s)",i,NaN);
	setResult("Flag",i,"not analyzed");
	updateResults();
	}

} else if (correction.length < 1) {
	roiManager("Select",i);
	roiManager("Rename","ROI_noPeak_"+i+1);
	print("ROI_"+i+1+"\r was not analyzed due to the difficulty to find the max fluorescence peak -high background");

	setResult("user Slices",i,NaN);
	setResult("corr Slices",i,NaN);
	setResult("Size-FWHM",i,NaN);
	setResult("decay -raw mean lifetime (frames)",i,NaN);
	setResult("decay -mean lifetime (s)",i,NaN);
	setResult("increase -raw mean lifetime(frames)",i,NaN);
	setResult("increase -mean lifetime (s)",i,NaN);
	setResult("Flag",i,"not analyzed");
	updateResults();
}

} else {
	roiManager("Select",i);
	print("ROI_"+i+1+"\r was not analyzed due to the bad fit for determining the event size -FWHM");

	setResult("user Slices",i,NaN);
	setResult("corr Slices",i,NaN);
	setResult("Size-FWHM",i,NaN);
	setResult("decay -raw mean lifetime (frames)",i,NaN);
	setResult("decay -mean lifetime (s)",i,NaN);
	setResult("increase -raw mean lifetime(frames)",i,NaN);
	setResult("increase -mean lifetime (s)",i,NaN);
	setResult("Flag",i,"not analyzed");
	updateResults();
	
}

} setBatchMode(false);

print("----------------");
//---------------------------------------------------------------------------------
// ----- Data extraction -----
number_ROI = roiManager("Count");
selectImage(ImageID);
run("Select None");
run("Duplicate...", "duplicate");
ImageID_final = getImageID();
selectImage(ImageID_final);
run("RGB Color"); 

for (n=0; n<number_ROI; n++) {
	size = getResult("Size-FWHM",n);
	decrease = getResult("decay -raw mean lifetime (frames)", n); 
	increase = getResult("increase -raw mean lifetime(frames)",n); 
	//increase = round(getResult("increase -raw mean lifetime(frames)",n)); 
	response = getResultString("Flag", n);
	max_fluorescence = getResult("Max Fluorescence",n);
	threshold_fluorescence = getResult("Fluorescence threshold",n);
	
	if (response == "not analyzed") {
		roiManager("Select",n); roiManager("Rename", "ROI_not analyzed_negative_"+n+1); setForegroundColor(255, 255, 0);
		getSelectionCoordinates(x_pos,y_pos);
		makeOval(x_pos[0]-size_condition/pixel_size, y_pos[0]-size_condition/pixel_size, 2*size_condition/pixel_size, 2*size_condition/pixel_size);
		if (Display ==true) {
			run("To Bounding Box");
			run("Enlarge...", "enlarge=" +manual_threshold*FWHM);
			
			b = getSliceNumber() - 4*increase_condition;
			if ( b < 1 ) { b = 1; }
			e = getSliceNumber() + 4*decrease_condition+1;
			if ( e > nSlices() ) { e = nSlices(); }
			run("Make Substack...", " slices="+b+"-"+e);
			run("8-bit"); rename("Event_ROI_not analyzed_negative"+n+1);
			saveAs("tiff", save_folder + name_file +"_Event_ROI_not analyzed_negative_"+n+1);
			close();
		} 
		selectImage(ImageID_final); run("Select None");
		makeOval(x_pos[0]-size_condition/pixel_size, y_pos[0]-size_condition/pixel_size, 2*size_condition/pixel_size, 2*size_condition/pixel_size);
		run("Draw", "slice");
		run("Select None");
		print("ROI_"+n+1);
		print("Event was not analyzed. See above for more details");
		
		

	} else {
		if (size >= size_condition && increase <= increase_condition && increase >= 0 && decrease < decrease_condition) {
			if (max_fluorescence >= threshold_fluorescence){
				roiManager("Select",n); roiManager("Rename", "ROI_positive_rescued_"+n+1); setForegroundColor(0, 255, 0);
				getSelectionCoordinates(x_pos,y_pos);
				makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
				if (Display ==true) {
					run("To Bounding Box");
					run("Enlarge...", "enlarge=" +manual_threshold*FWHM);
					
					b = getSliceNumber() - 4*increase_condition;
					if ( b < 1 ) { b = 1; }
					e = getSliceNumber() + 4*decrease_condition+1;
					if ( e > nSlices() ) { e = nSlices(); }
					run("Make Substack...", " slices="+b+"-"+e);
					run("8-bit"); rename("Event_ROI_positive_rescued_"+n+1);
					saveAs("tiff", save_folder + name_file +"_Event_ROI_positive_rescued_"+n+1);
					close();
				} 
				selectImage(ImageID_final); run("Select None");
				makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);

				run("Plot Z-axis Profile");
				Plot.getValues(x, y); close();
				str = "";
				for (i = 0; i < x.length ; i++) {
					str +=  "" + x[i] + "\t" + y[i] + "\n";
				}
				File.saveString( str, save_folder+"Intensity_Profile_positive_ROI_" + n+1 + ".xls" );
				
				selectImage(ImageID_final);
				run("Draw", "slice");
				run("Select None");
				print("ROI_"+n+1);
				print("Event was rescued due to a fluorescence intensity peak higher than Mean Fluorescence +" + "\r" + I_threshold + "xSD");
				
				
			} else {
				roiManager("Select",n); roiManager("Rename", "ROI_negative_"+n+1); setForegroundColor(255, 0, 0);
				getSelectionCoordinates(x_pos,y_pos);
				makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
				if (Display ==true) {
					run("To Bounding Box");
					run("Enlarge...", "enlarge=" +manual_threshold*FWHM); 
					
					b = getSliceNumber() - 4*increase_condition;
					if ( b < 1 ) { b = 1; }
					e = getSliceNumber() + 4*decrease_condition+1;
					if ( e > nSlices() ) { e = nSlices(); }
					run("Make Substack...", " slices="+b+"-"+e);
					run("8-bit"); rename("Event_ROI_negative_"+n+1);
					saveAs("tiff", save_folder + name_file +"_Event_ROI_negative_"+n+1);
					close();
				} selectImage(ImageID_final); run("Select None");
				makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
				run("Draw", "slice");
				run("Select None");
				print("ROI_"+n+1);
				
				print("Event does not fulfill the decrease condition and could not be rescued by fluorescence intensity");
				
			}
			if (response == "moving") {
				print("Event was not considered because it moved out of the confined boundary");
							
			}
			if (size < size_condition) {
				print("Event does not fulfill the size condition");
				
			}
			if (increase > increase_condition) {
				print("Event does not fulfill the fluorescence intensity increase condition");
				
			}
		} else {
			if (size >= size_condition && increase <= increase_condition && increase >= 0 && decrease <= number_slices && decrease >= decrease_condition) {
				if (response == "still") {
					roiManager("Select",n); roiManager("Rename", "ROI_positive_"+n+1); setForegroundColor(0, 255, 0);
					getSelectionCoordinates(x_pos,y_pos);
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
					if (Display ==true) {
						run("To Bounding Box");
						run("Enlarge...", "enlarge=" +manual_threshold*FWHM); 
						
						b = getSliceNumber() - 4*increase_condition;
						if ( b < 1 ) { b = 1; }
						e = getSliceNumber() + 4*decrease_condition+1;
						if ( e > nSlices() ) { e = nSlices(); }
						run("Make Substack...", " slices="+b+"-"+e);
						run("8-bit"); rename("Event_ROI_positive_"+n+1);
						saveAs("tiff", save_folder + name_file +"_Event_ROI_positive_"+n+1);
						close();
					} selectImage(ImageID_final); run("Select None");
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
					run("Plot Z-axis Profile");
					Plot.getValues(x, y);
					str = "";
					for (i = 0; i < x.length ; i++) {
						str +=  "" + x[i] + "\t" + y[i] + "\n";
					}
					File.saveString( str, save_folder+"Intensity_Profile_positive_ROI_" + n+1 + ".xls" );
					close(); selectImage(ImageID_final);
					
					run("Draw", "slice");
					run("Select None");
					print("ROI_"+n+1);
					print("Event fulfilled all the condition -no need for fluorescence rescue");
					
					
				}
				if (response == "moving") {
					roiManager("Select",n); roiManager("Rename", "ROI_negative_"+n+1); setForegroundColor(255, 0, 0);
					getSelectionCoordinates(x_pos,y_pos);
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
					if (Display ==true) {
						run("To Bounding Box");
						run("Enlarge...", "enlarge=" +manual_threshold*FWHM);
						
						b = getSliceNumber() - 4*increase_condition;
						if ( b < 1 ) { b = 1; }
						e = getSliceNumber() + 4*decrease_condition+1;
						if ( e > nSlices() ) { e = nSlices(); }
						run("Make Substack...", " slices="+b+"-"+e);
						run("8-bit"); rename("Event_ROI_negative_"+n+1);
						saveAs("tiff", save_folder + name_file +"_Event_ROI_negative_"+n+1);
						close();
					} selectImage(ImageID_final); run("Select None");
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);					
					run("Draw", "slice");
					run("Select None");
					print("ROI_"+n+1);
					print("Event was not considered because it moved out of the confined boundary");
					
									
				}
			} else {
				if (size < size_condition) {
					roiManager("Select",n); roiManager("Rename", "ROI_negative_"+n+1); setForegroundColor(255, 0, 0);
					getSelectionCoordinates(x_pos,y_pos);
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
					if (Display ==true) {
						run("To Bounding Box");
						run("Enlarge...", "enlarge=" +manual_threshold*FWHM); 
						
						b = getSliceNumber() - 4*increase_condition;
						if ( b < 1 ) { b = 1; }
						e = getSliceNumber() + 4*decrease_condition+1;
						if ( e > nSlices() ) { e = nSlices(); }
						run("Make Substack...", " slices="+b+"-"+e);
						run("8-bit"); rename("Event_ROI_negative_"+n+1);
						saveAs("tiff", save_folder + name_file +"_Event_ROI_negative_"+n+1);
						close();
					} selectImage(ImageID_final); run("Select None");	 	
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);			
					run("Draw", "slice");
					run("Select None");
					print("ROI_"+n+1); 
					print("Event does not fulfill the size condition");
					
						
				}
				if (increase > increase_condition) {
					roiManager("Select",n); roiManager("Rename", "ROI_negative_"+n+1); setForegroundColor(255, 0, 0);
					getSelectionCoordinates(x_pos,y_pos);
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
					if (Display ==true) {
						run("To Bounding Box");
						run("Enlarge...", "enlarge=" +manual_threshold*FWHM);
						
						b = getSliceNumber() - 4*increase_condition;
						if ( b < 1 ) { b = 1; }
						e = getSliceNumber() + 4*decrease_condition+1;
						if ( e > nSlices() ) { e = nSlices(); }
						run("Make Substack...", " slices="+b+"-"+e);
						run("8-bit"); rename("Event_ROI_negative_"+n+1);
						saveAs("tiff", save_folder + name_file +"_Event_ROI_negative_"+n+1);
						close();
					} selectImage(ImageID_final); run("Select None");
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);	
					run("Draw", "slice");
					run("Select None");
					print("ROI_"+n+1);
					print("Event does not fulfill the fluorescence intensity increase condition");
					
					 
				}
				if (increase < 0) {
					roiManager("Select",n); roiManager("Rename", "ROI_negative_"+n+1); setForegroundColor(255, 0, 0);
					getSelectionCoordinates(x_pos,y_pos);
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
					if (Display ==true) {
						run("To Bounding Box");
						run("Enlarge...", "enlarge=" +manual_threshold*FWHM);
						
						b = getSliceNumber() - 4*increase_condition;
						if ( b < 1 ) { b = 1; }
						e = getSliceNumber() + 4*decrease_condition+1;
						if ( e > nSlices() ) { e = nSlices(); }
						run("Make Substack...", " slices="+b+"-"+e);
						run("8-bit"); rename("Event_ROI_negative_"+n+1);
						saveAs("tiff", save_folder + name_file +"_Event_ROI_negative_"+n+1);
						close();
					} selectImage(ImageID_final); run("Select None");
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);	
					run("Draw", "slice");
					run("Select None");
					print("ROI_"+n+1);
					print("Event does not fulfill the fluorescence intensity increase condition: the increase is negative!");
					
						 
				}
				if (decrease > number_slices) {
					roiManager("Select",n); roiManager("Rename", "ROI_negative_"+n+1); setForegroundColor(255, 0, 0);
					getSelectionCoordinates(x_pos,y_pos);
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);
					if (Display ==true) {
						run("To Bounding Box");
						run("Enlarge...", "enlarge=" +manual_threshold*FWHM);
						
						b = getSliceNumber() - 4*increase_condition;
						if ( b < 1 ) { b = 1; }
						e = getSliceNumber() + 4*decrease_condition+1;
						if ( e > nSlices() ) { e = nSlices(); }
						run("Make Substack...", " slices="+b+"-"+e);
						run("8-bit"); rename("Event_ROI_negative_"+n+1);
						saveAs("tiff", save_folder + name_file +"_Event_ROI_negative_"+n+1);
						close();
					} selectImage(ImageID_final); run("Select None");
					makeOval(x_pos[0]-size/pixel_size, y_pos[0]-size/pixel_size, 2*size/pixel_size, 2*size/pixel_size);	
					run("Draw", "slice");
					run("Select None");
					print("ROI_"+n+1);
					print("Event does not fulfill the fluorescence intensity decrease condition: the decrease takes more time than the number of frames!");
					
			}
		}
	}
}
}	
print("----------------");
count_positive = 0; count_negative = 0; count_not_analyzed = 0;
number_ROI = roiManager("count");
for (i=0; i<number_ROI; i++) {
	roiManager("select",i); 
	label_count = getInfo("roi.name");
	if (label_count == "ROI_positive_"+i+1){
		count_positive +=1;
	}
	if (label_count == "ROI_positive_rescued_"+i+1){
		count_positive +=1;
	}
	if (label_count == "ROI_negative_"+i+1){
		count_negative +=1;
	}
	if (label_count == "ROI_not analyzed_negative_"+i+1){
		count_negative +=1;
		count_not_analyzed +=1;
	}
}

print("Total number of events:" +number_ROI);
print("Total number of positive events:" +count_positive);
print("Total number of negative events:" +count_negative);
print("(Total number of non analyzed events:" +count_not_analyzed+" counted as negative events)");

// Saving files
run("Select None");
roiManager("save", save_folder + name_file + "_Exocytosis Event XY location" + ".zip");

selectWindow("Log");
saveAs("Text", save_folder + name_file +"_Summary");
selectImage(ImageID_final);
saveAs("tiff", save_folder + name_file +"_Exocytosis Events");

//Save all the vesicles
selectWindow("Results");

lineseparator = "\n";
cellseparator = ",\t";

lines=split("Results", lineseparator);
labels=split(lines[0], cellseparator);
if (lengthOf(labels)==21) {
	k=5;
	for (j=k; j<labels.length; j++) {
		setResult(labels[j],0,0);
	}
	run("Clear Results");
	for (i=1; i<lines.length; i++) {
		parameters=split(lines[i], cellseparator);
        	for (j=k; j<parameters.length; j++){
        		setResult(labels[j],i-1,items[j]);
        		updateResults();
        	}
    }
}
saveAs("results", save_folder + name_file +"_Result table.csv");

// Save the positive vesicles
count_F = nResults-1;
for (j=count_F; j >= 0; j=j-1){
	if (getResultString("Flag", j) == "moving" || getResultString("Flag", j) == "not analyzed" || getResult("Size-FWHM",j) < size_condition || getResult("increase -raw mean lifetime(frames)",j) > increase_condition || getResult("increase -raw mean lifetime(frames)",j) < 0) {
	IJ.deleteRows(j, j);
	}
}
count_F = nResults-1;
for (j=count_F; j >= 0; j=j-1){
	if (getResult("decay -raw mean lifetime (frames)",j) < decrease_condition && getResult("Max Fluorescence",j) < getResult("Fluorescence threshold",j)) {	IJ.deleteRows(j, j);
	}
}
run("Summarize");
saveAs("results", save_folder + name_file +"_Result table_positive.csv");

// Cleaning
run("Close All");
open(save_folder + name_file +"_Exocytosis Events.tif");
setTool("point");

}





// ----------- FUNCTIONS -------------------------

// -- Extracting the rotation angle if R2 < 0.95 for the first try. This is to avoid a bad fit due to a huge burst
function IndexMax(array, value) {
	count=0; 
    for (i=0; i<lengthOf(array); i++) { 
        if (array[i]==value) { 
            count++; 
        } 
    }
    if (count>0) { 
        indices=newArray(count); 
        count=0; 
        for (i=0; i<lengthOf(array); i++) { 
            if (array[i]==value) { 
                indices[count]=i; 
                count++; 
            } 
        } 
        return indices;
    }
}

// -- Max value in an array
function maxOfArray(array) {
    min=0;
    for (c=0; c<lengthOf(array); c++) {
        min=minOf(array[c], min);
    }
    max=min;
    for (c=0; c<lengthOf(array); c++) {
        max=maxOf(array[c], max);
    }
    return max;
}

// -- Fitting outcome
function FitOutcome (x,y, curve, outcome) {
	Fit.doFit(curve, x, y);
	std_deviation = Fit.p(3);
	R2 = Fit.rSquared;
	FWHM = 2*pow(2*log(2),0.5)*std_deviation;
	if (outcome == "FWHM") {
		return FWHM;
	}
	if (outcome == "R2") {
		return R2;
	}
}