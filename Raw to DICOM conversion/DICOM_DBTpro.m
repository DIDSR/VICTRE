
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%               VICTRE (Virtual Imaging Clinical Trial for Regulatory Evaluation)

%          This matlab script converts DBT projection raw images (32-bit) from VICTRE to DICOM format.
 

%          Author: Eshan Dahal
%          eshan.dahal@fda.hhs.gov 

%          GITHUB LINK: https://github.com/DIDSR/VICTRE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

rawFolder = 'C:\raw-folder-path'; %specify the folder path with DBT projection images from VICTRE

dicom_folder = 'C:\dicom-folder-path'; %specify the folder path you want to write the DICOM images

filepattern = fullfile(rawFolder, 'mcgpu_image_pc_*.raw'); %example for signal absent data

Files = dir(filepattern);

for j = 1:length(Files)
    
    baseName = Files(j).name;
    
    num_ext1 = regexp(baseName, '\d*','match'); %extract patient no and projection number
    
    pro1 = cell2mat(num_ext1(2));

    pro_no1 = str2num(pro1); 
    
    if pro_no1 == 0
     
        index(j) = j;
        
        patient_id(j) = str2num(cell2mat(num_ext1(1))); %save unique patient ID   
    end  
end 
     

 index = index(index ~=0); % where 000 files are present 

 patient_id = patient_id(patient_id ~=0); %unique patients number 


study_uid= {1,length(index)};

series_uid =  {1,length(index)};

%remove 00 (FFDM) from Files for DBT dicom conversion 

for i = 1:length(index)
    
    study_uid{i} = {dicomuid}; % for each unique patient 
    
    series_uid{i}={dicomuid}; 
    
    u = index(i);

if i==1
    
    Files(u)=[ ];
    
else
    
    u = u-i+1;
    Files(u)=[ ];
    
end

end 


Files_2= Files; %truncated file 
    
    
Patient_ID = zeros(1,length(Files_2));
    

for kk=1:length(Files_2) 
    
    baseName1 = Files_2(kk).name;
     
    num_ext = regexp(baseName1, '\d*','match'); 
    
    Patient_ID(kk) = str2num(cell2mat(num_ext(1)));   
end

count = 0;

Patient_ID(max(kk)+1) = Patient_ID(max(kk)); 

for kkk=1:length(Files_2)    
     
    UID_test{kkk}= study_uid{1+count}; 
    
    SUID_test{kkk} = series_uid{1+count}; 
    
if  Patient_ID(kkk)~= Patient_ID(kkk+1) 
    
    UID_test{kkk+1} = study_uid{count+2}; 
    
    SUID_test{kkk+1} = series_uid{count+2}; 
    
     count = 1+ count;
end 

end 


for k=1:length(Files_2) 
    
    baseName1 = Files_2(k).name;
     
    num_ext = regexp(baseName1, '\d*','match');
     
    [a,b,c] = fileparts(baseName1);

    pro = cell2mat(num_ext(2));
 
    pro_no = str2num(pro); 
  
    fullName1 = fullfile(rawFolder, baseName1);
    
    fprintf(1, 'Now reading %s\n', fullName1);
    
    fid3 = fopen(fullName1, 'r');
    
    pCT = fread(fid3,3000*1500, 'float32'); 

    pmax = max(pCT); 
    
    pmin = min(pCT); 
   
    rescale = 8*(50+ (pCT-5200)*0.001469); %rescaling factor 
         
    cpCT = uint16(rescale);
     
    ipCT = reshape(cpCT,[3000 1500]);

    dicomwrite(ipCT, 'test.dcm');

    x = dicomread('test.dcm');

%%%%%%%%%%%%%%%%%META DATA DICOM TAG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%SOP CLass UID ****

metadata.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.1.2.1';

metadata.SOPClassUID = '1.2.840.10008.5.1.4.1.1.1.2.1';

%%%% Patient info %%%%%
metadata.PatientID = cell2mat(num_ext(1));  %patient ID 
metadata.PatientName = cell2mat(num_ext(1)); %same as patient ID
metadata.PatientComments = 'Dense breast type'; %Heterogenous, Scattered, Fatty
metadata.PatientState = 'Signal absent'; %Signal present or absent 
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
 metadata.PixelIntensityRelationship = 'LIN' ;
 metadata.PixelIntensityRelationshipSign = +1 ;  
 metadata.PresentationIntentType = 'FOR PROCESSING';
 metadata.PresentationLUTShape = 'IDENTITY';
 metadata.ReferringPhysicianName = ' ';
 metadata.RescaleIntercept = 0;
 metadata.RescaleSlope = 1;
 metadata.RescaleType = 'US';
 metadata.SeriesNumber = [];
 metadata.StudyDate = '20180106';
 metadata.StudyID = ' ';
 metadata.StudyTime = '1200';
 metadata.ViewCodeSequence = '';

% study section     
metadata.StudyDescription = 'Simulated DBT Projection';   
metadata.Modality = 'MG'; 
metadata.StudyInstanceUID = cell2mat(UID_test{k});

%%%%% SERIES under each modality %%%%%
metadata.SeriesDescription =  'Projection images';
metadata.SeriesInstanceUID = cell2mat(SUID_test{k});
metadata.AcquisitionNumber = pro_no; 
metadata.InstanceNumber = pro_no;
metadata.BodyPartExamined = 'BREAST'; 

%%%% Simulated Equipment %%%%%
 metadata.InstitutionName = 'FDA';
 metadata.InstitutionalDepartmentName = 'DIDSR';
 metadata.SoftwareVersions = 'MC-GPU_1.5b';
 
% %%%% Image %%%%%
 metadata.ImageType = 'ORIGINAL\PRIMARY';  
 metadata.ImagesInAcquisition = 25; 
 metadata.ImageComments = '85 x 85 micron pixel size; float32 to uint16 bit conversion'; 
 metadata.LossyImageCompression = '00'; 
 metadata.ConversionType = 'SYN'; 

% %%% Detector type %%%%%
 metadata.DetectorType = 'DIRECT' ; 
 metadata.DetectorConfiguration = 'AREA' ; 
 metadata.DetectorDescription = 'a-Se, 200 micron'; 
 metadata.DetectorActiveShape = 'RECTANGLE'; 
% 
% %%% X-ray Acquisition Dose %%%
 metadata.KVP = '28'; %depends on breast type
 metadata.ExposureInmAs = 3.5; 
 metadata.AnodeTargetMaterial = 'TUNGSTEN' ; 
 metadata.FilterType = 'FLAT' ; 
 metadata.FilterMaterial = 'RHODIUM'; 
 metadata.FilterThicknessMinimum = '0.050'; 
 
%%% Mammography image %%%%
metadata.DistanceSourceToDetector = '650'; % in mm from source to detector center 
metadata.DistanceSourceToPatient = '630' ; % in mm from source to the breast support side
metadata.PositionerType = 'MAMMOGRAPHIC'; 
metadata.DerivationDescription = 'rescaled raw data = 8*(50+ (pixel-5200)*0.001469)';
 
%%%%%%%%%%%%%%%%%%%%%%%% Write the metadata %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,b,c] = fileparts(baseName1);
name= '.dcm';
new_name =[num2str(b),name];
dicomwrite(x, fullfile(dicom_folder, new_name),metadata, 'CreateMode', 'copy');
fclose('all');

end 


