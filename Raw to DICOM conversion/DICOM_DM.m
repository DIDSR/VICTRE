

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%               VICTRE (Virtual Imaging Clinical Trial for Regulatory Evaluation)

%       This matlab script converts simulated digitial mamography raw images (32-bit) from VICTRE to DICOM format.%%% 

%          Author: Eshan Dahal
%          eshan.dahal@fda.hhs.gov 

%          GITHUB LINK: https://github.com/DIDSR/VICTRE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       							DISCLAIMER
%
%     This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. 
%     Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. 
%     Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, 
%     including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, 
%     and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, 
%     its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. 
%     Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. 
%     Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, 
%     and any modified versions bear some notice that they have been modified. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear

clc

close all

rawFolder ='C:\raw-folder-path'; %specify the folder path with digitial mamography images from VICTRE

dicom_folder = 'C:\dicom-folder';%specify the folder path you want to write the DICOM images

filepattern = fullfile(rawFolder, 'mcgpu_image_pc_*_0000.raw'); % example for signal absent data

Files = dir(filepattern);
  
for k=1:length(Files)
    
    baseName = Files(k).name;
    
    num_ext = regexp(baseName, '\d*','match'); %extract patient no and pro from the file name
    
    [a,b,c] = fileparts(baseName);
    
    fullName = fullfile(rawFolder, baseName);
    
    fprintf(1, 'Now reading %s\n', fullName);
    
    fid3 = fopen(fullName, 'r');
    
    pCT = fread(fid3,3000*1500, 'float32'); %reading DM VICTRE images

    pmax = max(pCT); 
    
    pmin = min(pCT); 
    
    rescale = 8*(50+ (pCT-5200)*0.000239);  %scaling factor     
    
    cpCT = uint16(rescale); 

    ipCT = reshape(cpCT,[3000 1500]);

    dicomwrite(ipCT, 'test1.dcm' );

    x = dicomread('test1.dcm');

%%%%%%%%%%%%%%%%%META DATA DICOM TAG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%SOP CLass UID ****

 metadata.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.1.2.1';

 metadata.SOPClassUID = '1.2.840.10008.5.1.4.1.1.1.2.1';

%%%% Patient info %%%%%
 metadata.PatientID = cell2mat(num_ext(1));  %patient ID 
 metadata.PatientName = cell2mat(num_ext(1)); %using same as patient ID
 metadata.PatientComments = 'Breast type info';  %Heterogenous, Scattered, Fatty
 metadata.PatientState = 'Signal absent or present'; %Signal present %Signal absent 
 metadata.PatientSize = []; 
 metadata.BodyPartThickness = 34.9; %in mm 


%%Clinical Trial Subject%%%
 metadata.ClinicalTrialProtocolName = 'VICTRE';
 metadata.ClinicalTrialSiteName = 'FDA';
	
%%%%Other information%%%%%%
 metadata.AccessionNumber = ' ';
 metadata.AcquisitionContextSequence = '';
 metadata.AnatomicRegionSequence = '';
 metadata.BurnedInAnnotation = 'NO';
 metadata.ClinicalTrialProtocolID = ' ';
 metadata.ClinicalTrialSiteID = ' '; 
 metadata.ClinicalTrialSponsorName = ' ';
 metadata.ClinicalTrialSubjectID = ' ';
 metadata.ClinicalTrialSubjectReadingID = ' ';
 metadata.ImageLaterality= 'R';
 metadata.ImagerPixelSpacing = '0.085\0.085';
 metadata.InstanceNumber = '';
 metadata.Manufacturer = 'VICTRE';
 metadata.OrganExposed = 'BREAST';
 metadata.PatientBirthDate = '';
 metadata.PatientOrientation = 'P\H';
 metadata.PatientSex = 'F';
 metadata.PixelIntensityRelationship = 'LIN' ; %LIN or LOG 
 metadata.PixelIntensityRelationshipSign = +1 ;  
 metadata.PresentationIntentType = 'FOR PROCESSING';
 metadata.PresentationLUTShape = 'IDENTITY';
 metadata.ReferringPhysicianName = ' ';
 metadata.RescaleIntercept =  0;
 metadata.RescaleSlope = 1;
 metadata.RescaleType = 'US';
 metadata.SeriesNumber = [];
 metadata.StudyDate = '20180106';
 metadata.StudyID = ' ';
 metadata.StudyTime = '1200';
 metadata.ViewCodeSequence = '';

%%%% STUDY in terms of Modalities %%%%% 

 metadata.StudyDescription = 'Simulated Digital Mammography';  
 metadata.Modality = 'MG'; 
 metadata.StudyInstanceUID = dicomuid; %Study instance UID automatically generated

%%%%% SERIES under each modality %%%%%

 metadata.SeriesDescription = num2str(b);  % File name of the series                                                                                                                                                              
 metadata.BodyPartExamined = 'BREAST'; 
 metadata.SeriesInstanceUID = dicomuid;  %Series instance UID automatically generated                                                                                                                                                     metadata.BodyPartExamined = 'BREAST'; 

%%%% Simulated Equipment %%%%%
 metadata.InstitutionName = 'FDA';
 metadata.InstitutionalDepartmentName = 'DIDSR';
 metadata.SoftwareVersions = 'MC-GPU_1.5b';
 
% %%%% Image %%%%%
 
 metadata.ImageType = 'ORIGINAL\PRIMARY'; 
 metadata.ImageComments = '85 x 85 micron pixel size; 5:1 antiscatter grid used; float 32 to uint16 bit conversion'; 
  
% %metadata.AcquisitionNumber = ' '; 
 metadata.ImagesInAcquisition = 1; 
 metadata.LossyImageCompression = '00';
 metadata.ConversionType = 'SYN'; %


% %%% Detector type %%%%%
 metadata.DetectorType = 'DIRECT' ; 
 metadata.DetectorConfiguration = 'AREA' ; 
 metadata.DetectorDescription = 'a-Se, 200 micron'; 
 metadata.DetectorActiveShape = 'RECTANGLE'; 
% 
% %%% X-ray Acquisition Dose %%%
 metadata.KVP = '28'; 
 metadata.ExposureInmAs = 58.6;
 metadata.AnodeTargetMaterial = 'TUNGSTEN' ; 
 metadata.FilterType = 'FLAT' ; 
 metadata.FilterMaterial = 'RHODIUM'; 
 metadata.FilterThicknessMinimum = '0.050'; %in mm

%%% Mammography image %%%%
 metadata.DistanceSourceToDetector = '650'; % in mm from source to detector center 
 metadata.DistanceSourceToPatient = '630' ; % in mm from source to the breast support side
 metadata.PositionerType = 'MAMMOGRAPHIC';
 
 metadata.DerivationDescription = 'rescaled raw data = 8*(50+ (pixel-5200)*0.000239)';

%%%%%%%%%%%%%%%%%%%%%% write the metadata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [a,b,c] = fileparts(baseName);

 name= '.dcm';

 new_name =[num2str(b),name];

 dicomwrite(x, fullfile(dicom_folder, new_name),metadata, 'CreateMode', 'copy');

 fclose('all');

end 





