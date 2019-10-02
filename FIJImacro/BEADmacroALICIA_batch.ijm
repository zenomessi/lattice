print("\\Clear"); //clear log window
print("====== BEADmacro ======"); //title
//
requires("1.45s"); //according to my ImageJ ;)
//*******************************************************************
//Startup GUI (choose what to do: filter/background, align, PIV, Map)
//*******************************************************************
/*Dialog.create("BEADmacro");
Dialog.addCheckbox("(1) Pretreat",false);
Dialog.addCheckbox("(2) Duplicate",false);
Dialog.addCheckbox("(3) Realign",false);
Dialog.addCheckbox("(4) PIV",false);
Dialog.addCheckbox("(5) Map",false);
Dialog.show();
do_pretreat = Dialog.getCheckbox();
do_duplicate = Dialog.getCheckbox();
do_realign = Dialog.getCheckbox();
do_piv = Dialog.getCheckbox();
do_map = Dialog.getCheckbox();
if( !(do_pretreat|do_duplicate|do_realign|do_piv|do_map) ){ //nothing to do?
	print("Nothing to do! End.");
	exit;
}*/
//********
//Open GUI
//********
#@String wdir
#@String filename
print("=> Choose working directory..."); //log
//currentDirectory = getDirectory("Choose working directory"); //choose dir
currentDirectory = wdir
print("Working directory: "+wdir);  //log
print("=> Choose a file..."); //log
/*fileList = getFileList(currentDirectory); //list available files
Dialog.create("Pick a .tif file:"); //show available files
Dialog.addChoice("File",fileList); //comboBox
Dialog.show();
fileName = Dialog.getChoice(); //retrieve choice*/
print("File: "+filename); //log
print(currentDirectory+filename);
if( File.isDirectory(currentDirectory+filename) ){ //error: directory instead of file
	beep();
	print("ERROR: Choice is a directory! End."); //log
	exit("ERROR: Choice is a directory! End.");
}
if( !endsWith(filename,"tif") ){ //error: wrong file type
	beep();
	print("ERROR: Not a Tiff file! End."); //log
	exit("ERROR: Not a Tiff file! End.");
}
//**************
//AutoDimensions
//**************
open(currentDirectory+filename); //image open
if( nSlices <=1 ){ //Error check: less than 2 slices?
	beep();
	print("ERROR: Image is not a stack! End."); //log
	exit("ERROR: Image is not a stack! End.");
}
width = getWidth(); //image width
height = getHeight(); //image height
//REMARK: <nSlices> is automatically set when opening an image
print("Image is a "+width+" by "+height+" stack, "+nSlices+" slices."); //log
//***************
//Output settings
//***************
pathOut  = currentDirectory+"pivout"+File.separator;
print("Output dir: "+pathOut); //log
if( File.exists(pathOut) ) { //output directory already exists
	print("Using existing /pivout directory...");
} else { //output directory does NOT exist
	print("Creating /pivout subdir...");
	File.makeDirectory(pathOut); //make it
	if( !File.exists(pathOut) ){ //error making dir
		beep();
		print("ERROR: Unable to create output directory! End.");
		exit("ERROR: Unable to create output directory! End.");
	}
}
//************
//(1) PRETREAT
//************
//if( do_pretreat == true ){
	print("PRETREAT..."); //log
	//pretreat GUI
	//------------
	//print("=> Choose pretreatment parameters..."); //log
	print("=> Pretreatment parameters..."); //log
	print("Median Filter Radius 3")
	print("Length of sliding window 20")
	/*Dialog.create("Filter and Background");
	Dialog.addCheckbox("Median filter",true);
	Dialog.addNumber("Radius",3);
	Dialog.addCheckbox("Background subtract",true);
	Dialog.addNumber("Length of sliding window",20);
	Dialog.show();
	do_median = Dialog.getCheckbox();
	median_radius = Dialog.getNumber;
	do_background = Dialog.getCheckbox();
	background_slidingWindow = Dialog.getNumber;
	//median filter
	//-------------
	if( do_median == true ){ //median filter selected
		print("Median filtering (radius = "+median_radius+" px)...");
		run("Median...", "radius="+median_radius+" stack");
	}*/
	run("Median...", "radius=3 stack");
	run("Background Subtractor", "length=20 stack");
	//background subtraction
	//----------------------
	if( do_background == true ){ //background subtract selected
		if( File.exists(getDirectory("plugins")+"mosaic_plugins.jar") ){ //try to use MOSAIC plugin
			print("Subtracting background with MOSAIC plugin (length of sliding window ="+background_slidingWindow+"px)...");
			run("Background Subtractor", "length="+background_slidingWindow+" stack");
		}
		else{ //otherwise rollback to ij default background subtract
			print("Subtracting background using ImageJ default (length of sliding window ="+background_slidingWindow+"px)...");
			run("Subtract Background...", "rolling="+background_slidingWindow+" sliding stack"); //sliding paraboloid + smoothing enabled
		}
	}
	if( do_median | do_background ){ //some pretreatment done?
		run("Enhance Contrast", "saturated=0.35"); //try some auto brightness/contrast
		resetMinAndMax();
		fileName = "t_"+fileName; //upgrade to use the treated image in the following, if any
		print("Saving to "+fileName+"..."); //log
		save(currentDirectory+fileName); //save image as t(reated)_<oldname>
		close();
		open(currentDirectory+fileName);

	}
//}

//************
//(2) DUPLICATE
//************
if( do_duplicate == true ){
	print("DUPLICATE..."); //log
	i=indexOf(fileName,".tif");
	fileNameNoExtension = substring(fileName,0,i);
	run("Duplicate...", "duplicate range="+nSlices+"-"+nSlices);
	run("Concatenate...", "image1="+fileNameNoExtension+"-1.tif"+" image2="+fileName+" image3=[-- None --]");
	setSlice(nSlices);
	run("Delete Slice");
	//setSlice(1);
	fileName = "d_"+fileName; //upgrade to use the treated image in the following, if any
	print("Saving to "+fileName+"..."); //log
	save(currentDirectory+fileName); //save image as d(uplicate)_<oldname>
	close();
	if( do_piv ){
		open(currentDirectory+fileName);
	}
}
//***********
//(3) REALIGN
//***********
if( do_realign==true ){
	print("REALIGN..."); //log
	if( !File.exists(getDirectory("plugins")+"Template_Matching.jar") ){ //plugin available?
		beep();
		print("ERROR: Missing plugin (<Template_Matching.jar>)! End.");
		exit("ERROR: Missing plugin (<Template_Matching.jar>)! End.");
	}
	//realign
	//-------
	print("=> Set alignment parameters and draw a ROI for pattern matching..."); //log
	run("Align slices in stack...", "");
	//++
	roiManager("Add");
	roiManager("Save",currentDirectory+"roi.zip");
	selectWindow("ROI Manager");
	run("Close");
	//++
	//evaluate displacement
	//---------------------
	displacement_array = newArray(nResults); //create an array for displacement vector =[(dX^2+dY^2)^0.5]
	for( i=0; i<nResults; i+=1 ){ //foreach slice...
		displacement_array[i] = sqrt( pow(getResult("dX",i),2) + pow(getResult("dY",i),2) ); //compute displacement vector
	}
	Array.getStatistics(displacement_array, min, max, mean, std); //compute min mean max displacement
	print("Average displacement vector = "+mean+"px (min "+min+"px, max "+max+"px)."); //log
	if( mean <3 ){ //average realignment smaller than 3px
		print("WARNING: average displacement vector <3px!"); //warn
	}
	selectWindow("Results");
	saveAs("measurements",currentDirectory+"realign_results.txt"); //save displacement table
	run("Close"); //and close it
	makeRectangle(0,0,0,0); //delete pattern matching roi
	fileName = "r_"+fileName; //upgrade to use the treated image in the following, if any
	print("Saving to "+fileName+"..."); //log
	save(currentDirectory+fileName); //save image as r(ealigned)_<oldname>
	close();
	if( do_duplicate ){
		open(currentDirectory+fileName);
	}
}



//************
//(4) PIV
//************
if( do_piv==true ){ //ex "MacroForPIV-advanced-Realigned.ijm"
	print("PIV...");
	if( !File.exists(getDirectory("plugins")+"PIV_.jar") ){
		beep();
		print("ERROR: Missing plugin (<PIV_.jar>)! End.");
		exit("ERROR: Missing plugin (<PIV_.jar>)! End.");
	}
	//parameter GUI (default suggested)
	Dialog.create("PIV parameters");
	Dialog.addNumber("PIV pass #1 window",250);
	Dialog.addNumber("PIV pass #1 search window",350);
	Dialog.addNumber("PIV pass #1 vector size",125);
	Dialog.addNumber("PIV pass #2 window",100);
	Dialog.addNumber("PIV pass #2 search window",250);
	Dialog.addNumber("PIV pass #2 vector size",50);
	Dialog.addNumber("PIV pass #3 window",50);
	Dialog.addNumber("PIV pass #3 search window",150);
	Dialog.addNumber("PIV pass #3 vector size",35);
	Dialog.addNumber("Correlation threshold",0.9);
	Dialog.addNumber("Pixel/um",0.1075);
	Dialog.addNumber("Young modulus (Pa)",1500);
	Dialog.addNumber("Poisson ratio",0.5);
	Dialog.addNumber("Regularization",0.0000000009);
	Dialog.show();
	piv1 = Dialog.getNumber();
	sw1 = Dialog.getNumber();
	vs1 = Dialog.getNumber();
	piv2 = Dialog.getNumber();
	sw2 = Dialog.getNumber();
	vs2 = Dialog.getNumber();
	piv3 = Dialog.getNumber();
	sw3 = Dialog.getNumber();
	vs3 = Dialog.getNumber();
	correlation = Dialog.getNumber();
	pixel_um = Dialog.getNumber();
	young = Dialog.getNumber();
	poisson = Dialog.getNumber();
	regularization = Dialog.getNumber();
	//go!
	run("8-bit");
	for (n=2; n<=nSlices; n++){
		run("Make Substack...", "  slices=1,"+n);
		if (n<10){nn="0"+n;} else{nn=n;}
		outn=pathOut+"PIV_Stack"+nn+".txt";
		run("iterative PIV(Advanced)...", "piv1="+piv1+" sw1="+sw1+" vs1="+vs1+" piv2="+piv2+" sw2="+sw2+" vs2="+vs2+" piv3="+piv3+" sw3="+sw3+" vs3="+vs3+" correlation="+correlation+" use debug_x=-1 debug_y=-1 path="+outn);
		wait( 100 );
		run("FTTC ", "pixel="+pixel_um+" poisson="+poisson+" young's="+young+" regularization="+regularization+" plot="+width+" plot="+height+" select="+outn);
//		run("FTTC ", "pixel=0.107 poisson=0.5 young's=2500 regularization=0.000000001 plot=1392 plot=1392 select="+outn);
		run("Close");
		run("Close");
		run("Close");
		run("Close");
		run("Close");
//		run("Close");
//		run("Close");
//		run("Close");
	}
	selectWindow("Results");
	run("Close"); //and close it
	if( !do_map ){
		close();
	}
}

//************
//(5) MAP
//************
if( do_map==true ){
	print("MAP...");
	list = getFileList(pathOut);
	if( list.length<1 ){
		exit("ERROR: no .txt file to process! End.");
	}
	setBatchMode(true);
	number_of_slices = nSlices;
	close();
	run("Hyperstack...", "title=ForceDisplacement type=16-bit display=Color width="+width+" height="+height+" channels=2 slices=2 frames="+(number_of_slices-1)+" label");
	print("Hyperstacking:");
	for (n=2; n<number_of_slices; n+=1){
		if (n<10){nn="0"+n;} else{nn=n;}
		print("Slice "+n+"...");
		dispFile=pathOut+"PIV_Stack"+nn+".txt";
		forceFile=pathOut+"Traction_PIV_Stack"+nn+".txt";
		run("plot...", "select="+dispFile+" select="+dispFile+" vector_scale=1 max=2 plot_width=0 plot_height=0 show draw lut=S_Pet");
		run("plot FTTC", "select="+forceFile+" select="+forceFile+" vector_scale=0.01 max=2000 plot_width=0 plot_height=0 show draw lut=S_Pet");
		selectWindow("Scale Graph");
		close();
		selectWindow("Scale Graph");
		close();
		selectWindow("Vector plot_PIV_Stack"+nn+".txt");
		run("Copy");
		close();
		selectWindow("ForceDisplacement");
		setSlice((n-1)*4+1);
		run("Paste");
		selectWindow("Magnitude map_PIV_Stack"+nn+".txt");
		run("Copy");
		close();
		selectWindow("ForceDisplacement");
		setSlice((n-1)*4+2);
		run("Paste");
		selectWindow("Vector plot_Traction_PIV_Stack"+nn+".txt");
		run("Copy");
		close();
		selectWindow("ForceDisplacement");
		setSlice((n-1)*4+3);
		run("Paste");
		selectWindow("Magnitude map_Traction_PIV_Stack"+nn+".txt");
		run("Copy");
		close();
		selectWindow("ForceDisplacement");
		setSlice((n-1)*4+4);
		run("Paste");
	}
	print("Hyperstacking Done!");
	setBatchMode("exit and display");
	Stack.setChannel(2);
	Stack.setSlice(1);
	run("S_Pet");
	save(currentDirectory+"h_"+fileName); //save image as h_(yperstack)_<oldname>
	//try brightness/contrast autoset...
	Stack.setChannel(2);
	Stack.setSlice(1);
	Stack.setFrame( round(number_of_slices/2) ); //... in the middle of the sequence
	run("Enhance Contrast", "saturated=0.35");
	Stack.setFrame(1); //go home
	Stack.setFrameRate(2);
	doCommand("Start Animation [\\]"); //and play
}
print("====== BEADmacro ==end="); //title
selectWindow("Log");
saveAs("text",currentDirectory+"log.txt"); //save log
