
file name: countnuclei.m

Procedures to run the file
1. Make sure to copy and paste the root cell images into the same folder of this file 

2. type the following in Matlab editor:
	img1 = imread('StackNinja1.bmp');
	result = countnuclei(img1);

3. you shall see different sub-images pop up in different figures in sequence and with descriptions.
	it would take awhile before you can see the plots of the nuclei size, shape and brightness distribution
	because of the for-loop
	
4. When the program is done running,in the Matlab editor, check out the variable "result" in the workspace.
	
5. inside the struct "result", there will be a list of related variables available for view.

	to view the final extracted nuclei in binary, click "result.combine_binary"
	to view the final extracted nuclei in color, click "result.combine_color"
	to view the total count of nuclei, click "result.sum_nuclei"
	
	to view the size data of each nucleus region, click  "result.nucleiarea"
	to view the shape data of each nucleus region, click "result.nucleisolidity"
														 "result.nucleieccentricity"
														 "result.nuclei_bound"
														 "result.nuclei_chaincode"
	
	to view the brightness of each nucleus region, click "result.nucleus_brightness"
	to view the sum of brightness of all nuclei regions, click "result.sum_nuclei_brightness"
	to view the average brightness of all nuclei regions, click "result.nuclei_avgbrightness"
	
6. Repeat step 1 to 5 for 
	img2 = imread('StackNinja2.bmp');
	img3 = imread('StackNinja3.bmp');
	

