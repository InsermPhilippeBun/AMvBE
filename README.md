# AMvBE
CD63- and CD81-pHluorin report multivesicular body (MVB)-plasma membrane (PM) fusion events that appear as sudden fluorescently increasing spots that decrease over time. AMvBE is an ImageJ/FIJI macro that provides a non-biased analysis of these events by following criteria initially determined by users. The macro comprises two steps: (1) a manual step that involves pinpointing events that resemble fusion events and (2) an automatic analysis of chosen events to extract parameters (decrease and increase mean lifetimes, fluorescent intensity, displacement and MVB estimated size) which are further tested against MVB-PM fusion requirement criteria. 

How to install and run macro
Download ImageJ or FIJI on your computer first. Ensure that additional plugins (MultiStackReg and TurboReg) are installed. If not, download and install the additional plugins (see below).
AMvBE macro has been successfully tested on ImageJ 1.50 version or higher on both MAC OSX and Windows 10 environments with JAVA 7 or higher. The typical installation time is less than 5 minutes.
Download the macro file “AMvBE.ijm” and store it on your hard drive.
To start the macro, install the macro (Plugins/Macro/Install…) and Press Run (or Ctrl + R)
Additional plugin requirements
The macro uses the MultiStackReg plugin (version 1.45 is available here - http://bradbusse.net/downloads.html) to register stacks of fluorescent images. This plugin requires the use of TurboReg, which can be found here (http://bigwww.epfl.ch/thevenaz/turboreg/ - Thevenaz et al., IEEE Transactions on Image Processing, 1998). Read the ImageJ documentation to learn how to install plugins. In short, place the downloaded MultiStackReg1.45_.jar and TurboReg_.jar files in the Plugins folder of FIJI. Then start FIJI and install TurboReg (Plugins/Install Plugin/select and open TurboReg_.jar), after successful installation restart FIJI. Repeat the same procedure for MultiStackReg and restart FIJI again. After restarting, the plugins should appear in the Plugin list (MultiStackReg should appear under Plugins/Registration/MultiStackReg). 
Bleach correction can be performed using the plugin already installed on ImageJ/FIJI (Image/Adjust/Bleach Correction). Details can be found here.
Analysis parameters
MVB-PM fusion events, and as a corollary exosome release, are visible as a burst of fluorescent intensty from the pHluorin tag. To count and sort these events, users can tune the recommended settings of the analysis parameters, which are detailed in Box 2 of Bebelman et al., Nat Prot, 2019 and below (see Bebelman et al., Nat Prot, 2019 Fig 3c for an explanatory schematic of the parameters).
-Time window for event correction (in nb of frames): 6
During the selection of the events, the user might pinpoint a fusion event some frames before or after the peak fluorescence intensity of the event is reached. The time window for event correction allows AMvBE to determine the maximum fluorescence intensity of the event 6 frames before and after the frame pinpointed by the user (Fig 3c, 1,2). The higher this value, the less accurate the user has to pinpoint the right frame. However, higher values might result in problems when there are multiple fusion events that occur close together in space and time. 
-Condition for fluorescence decrease after exocytosis (min. nb of frames): 3
The fluorescence decay after fusion is fitted as exponential decay with I(t)=Imax*e^-(t/τ), where I = intensity and τ is the decay rate. The duration of the fluorescence decay is the time from peak fluorescence intensity until the timepoint where  I(t)=Imax*e^-((t=τ)/τ)), which corresponds to the timepoint where the intensity is ±0.37 times the maximum intensity (Fig 3c, 3-5,8,9). 
-Condition for fluorescence increase before exocytosis (max. nb of frames): 3
The fluorescence increase before fusion is fitted as exponential growth with I(t)=Imax*e^(t/τ), where I = intensity and τ is the increase rate. The duration of the fluorescence increase is the time from the timepoint where I(t)=Imax*e^((t=τ)/τ)) until the peak fluorescent intensity is reached (Fig 3c, 3-7).
-Min. event size to be considered (in micrometer): 0.4
The size of the fusion event (vesicle diameter) is measured as the full width at half the maximum intensity of the fusion event (Fig 3c, 10).
-If the decay condition is not fulfilled, threshold value (N) for exocytosis event rescue such that Mean + N*STD: 3
The fusion of MVBs that contain a relatively high TSPAN-pHluorin concentration in the limiting membrane can result in a short fluorescence decrease time until   is reached, followed by a slower residual decrease of the remaining exosome-associated TSPAN-pHluorin fluorescence. Optionally, AMvBE can rescue these MVB-PM fusion events by setting this parameter at 3. If deemed inappropriate, this parameter can be set at 100 to avoid rescue of these events.
-Manual threshold for vesicle movement fixed by users (unit in nb of vesicle diameter): 2
To avoid counting kiss-and-run-like fusion events, the maximum movement of the fluorescent signal after fusion is restricted to twice the vesicle diameter (Fig 3c, 10,11)

Output
Upon analysis of fluorescence bursts, a text file named “Results” is generated. This document contains all the user-defined parameters as well as the event sorting outcome which comes as followed:  - positive events that fulfilled all user-defined criteria  - negative events that did not fulfill user-defined criteria. The failed criterion (-a) is (are) mentioned.  - positive rescued events that fulfilled all criteria except for the decrease condition but rescued by the fluorescence intensity peak condition.  - not analyzed events that failed to be properly analyzed. Further explanations are given. 
In addition, pinpointed locations of fluorescent bursts (supposedly MVB fusion events) before and after spatiotemporal correction are saved as well as Mean fluorescence intensity profile of each positive event and result table files (csv format). Finally, a RGB image is generated with colored circle with respect to the outcome of event sorting (red: negative; yellow: not analyzed events; green: positive events). Optionally users can decide to generate stacks centered on positive events.
