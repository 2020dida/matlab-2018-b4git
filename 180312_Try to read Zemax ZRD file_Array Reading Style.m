clear all
fclose('all')
%%
File_Path='C:\TuanShu\180312_MLOPTIC\180312_SL Test_All Lens 99_1000 rays_100W_Final deisgn of DET_Change Element order.ZRD';
% 1st int32: version
% 2nd int32: number of segment (max)
% 3rd in32: number of segment of ray 1
Known_Number_of_Ray=1000;
Estimated_Number_of_Segment=300;
Size_of_Segment=208;    %byte
Size_In32=Size_of_Segment/4;
fin=fopen(File_Path);
Version=fread(fin,1,'int32');
MaxNSeg=fread(fin,1,'int32');

RayNumber=1;
NumberofSegmentforThisRay=zeros(Known_Number_of_Ray,1);
RayData=cell(Estimated_Number_of_Segment,Known_Number_of_Ray);
%RayData=struct('status',{},'level',{},'hit_object',{},'hit_face',{},'unused',{},'in_object',{},'parent',{},'storage',{},'xybin',{},'lmbin',{},'index',{},'starting_phase',{},'x',{},'y',{},'z',{},'l',{},'m',{},'n',{},'nx',{},'ny',{},'nz',{},'path_to',{},'intensity',{},'phase_of',{},'phase_at',{},'exr',{},'exi',{},'eyr',{},'eyi',{},'ezr',{},'ezi',{});   %there are 31 type of properties of a segment
while RayNumber<=Known_Number_of_Ray
    NumberofSegmentforThisRay(RayNumber)=fread(fin,1,'int32');
    if NumberofSegmentforThisRay(RayNumber) == []
        RayNumber=0;
        break;
    end
    disp(RayNumber);
    for SegNumber=1:NumberofSegmentforThisRay(RayNumber)  % Segment RayNumber
        Temp1=fread(fin,10,'uint32');    %read 1~10
        

%         unsigned int status; 1
        RayData{SegNumber,RayNumber}.status=Temp1(1);

%         int level;        2
        RayData{SegNumber,RayNumber}.level=Temp1(2);

%         int hit_object;   3
        RayData{SegNumber,RayNumber}.hit_object=Temp1(3);

%         int hit_face; 4
        RayData{SegNumber,RayNumber}.hit_face=Temp1(4);

%         int unused;   5
        RayData{SegNumber,RayNumber}.unused=Temp1(5);

%         int in_object;    6
        RayData{SegNumber,RayNumber}.in_object=Temp1(6);

%         int parent;   7
        RayData{SegNumber,RayNumber}.parent=Temp1(7);

%         int storage;  8
        RayData{SegNumber,RayNumber}.storage=Temp1(8);

%         int xybin, lmbin; 9 10
        RayData{SegNumber,RayNumber}.xybin=Temp1(9);
        RayData{SegNumber,RayNumber}.lmbin=Temp1(10);

        Temp2=fread(fin,21,'double');    %read 11~31 (total 21)

%         double index, starting_phase; 11 12
        RayData{SegNumber,RayNumber}.index=Temp2(1);
        RayData{SegNumber,RayNumber}.starting_phase=Temp2(2);

%         double x, y, z;    13 14 15

        RayData{SegNumber,RayNumber}.x=Temp2(3);
        RayData{SegNumber,RayNumber}.y=Temp2(4);
        RayData{SegNumber,RayNumber}.z=Temp2(5);

%         double l, m, n;    16 17 18

        RayData{SegNumber,RayNumber}.l=Temp2(6);
        RayData{SegNumber,RayNumber}.m=Temp2(7);
        RayData{SegNumber,RayNumber}.n=Temp2(8);

%         double nx, ny, nz;     19 20 21

        RayData{SegNumber,RayNumber}.nx=Temp2(9);
        RayData{SegNumber,RayNumber}.ny=Temp2(10);
        RayData{SegNumber,RayNumber}.nz=Temp2(11);

%         double path_to, intensity;     22 23
        RayData{SegNumber,RayNumber}.path_to=Temp2(12);
        RayData{SegNumber,RayNumber}.intensity=Temp2(13);

%         double phase_of, phase_at;      24 25

        RayData{SegNumber,RayNumber}.phase_of=Temp2(14);
        RayData{SegNumber,RayNumber}.phase_at=Temp2(15);

%         double exr, exi, eyr, eyi, ezr, ezi;   26 27 28 29 30 31

        RayData{SegNumber,RayNumber}.exr=Temp2(16);
        RayData{SegNumber,RayNumber}.exi=Temp2(17);

        RayData{SegNumber,RayNumber}.eyr=Temp2(18);
        RayData{SegNumber,RayNumber}.eyi=Temp2(19);

        RayData{SegNumber,RayNumber}.ezr=Temp2(20);
        RayData{SegNumber,RayNumber}.ezi=Temp2(21);

        %fseek(fin,Size_of_Segment, 'cof');
        %disp(SegNumber);
    end
    RayNumber=RayNumber+1;
end

%%
SegNumber=189;
RayNumber=16;
RayData{SegNumber,RayNumber}.hit_object
%% Check the summation of the intensity of the last segment of every rays
Sum=0;
for RayNumber=1:Known_Number_of_Ray
    Sum=Sum+RayData{NumberofSegmentforThisRay(RayNumber),RayNumber}.intensity;
    disp(RayNumber);
end
%% Check the max (last) hit surface of each ray

LastSurface=zeros(Known_Number_of_Ray,1);
for RayNumber=1:Known_Number_of_Ray
    for SegNumber=1:NumberofSegmentforThisRay(RayNumber)
        if RayData{SegNumber,RayNumber}.hit_object>LastSurface(RayNumber);
            LastSurface(RayNumber)=RayData{SegNumber,RayNumber}.hit_object;
        end
    end
        disp(RayNumber);
end

%% Check the ray path of specific ray
RayNumber=392;
Surface_Array=zeros(NumberofSegmentforThisRay(RayNumber),1);
Level_Array=zeros(NumberofSegmentforThisRay(RayNumber),1);
Status_Array=zeros(NumberofSegmentforThisRay(RayNumber),1);
for SegNumber=1:NumberofSegmentforThisRay(RayNumber)
    Surface_Array(SegNumber)=RayData{SegNumber,RayNumber}.hit_object;
    Level_Array(SegNumber)=RayData{SegNumber,RayNumber}.level;
    Status_Array(SegNumber)=RayData{SegNumber,RayNumber}.status;

end

plot(1:NumberofSegmentforThisRay(RayNumber),Surface_Array,1:NumberofSegmentforThisRay(RayNumber),Level_Array,1:NumberofSegmentforThisRay(RayNumber),Status_Array);
%% Identify rays, trace parent, check the last (max) surface intersect
Det_Surface=15;
% 先找hit Detector plane之segment
% 再trace他的parant
Count=0;
Error_Count=0;
ROI_LastOBJ_Array=[];
%ROI_LastSurf_Array=[];
ROI_Intensity_at_DET_Array=[];
for RayNumber=1:Known_Number_of_Ray
    for SegNumber=1:NumberofSegmentforThisRay(RayNumber)
        if RayData{SegNumber,RayNumber}.hit_object==Det_Surface;
            Current_Seg=SegNumber;
            Current_Obj=Det_Surface;
            %Current_Surf=RayData{Current_Seg,RayNumber}.hit_face;
            Count=Count+1;
            ROI_Intensity_at_DET_Array(Count)=RayData{Current_Seg,RayNumber}.intensity;
            while Current_Seg
                Current_Seg=RayData{Current_Seg,RayNumber}.parent;
                if Current_Seg==0
                    break;
                end
                Current_Obj=RayData{Current_Seg,RayNumber}.hit_object;
                %Current_Surf=RayData{Current_Seg,RayNumber}.hit_face;
                if Current_Obj==Det_Surface
                    break;
                    Count=Count-1;
                    Error_Count=Error_Count+1;
                end
                if length(ROI_LastOBJ_Array)<Count
                    ROI_LastOBJ_Array(Count)=Current_Obj;
                %ROI_LastSurf_Array(Count)=Current_Surf;

                elseif Current_Obj>ROI_LastOBJ_Array(Count)
                    ROI_LastOBJ_Array(Count)=Current_Obj;
                    %ROI_LastSurf_Array(Count)=Current_Surf;

                end       
            end
        end
    end
        disp(RayNumber);
end
histogram(ROI_LastOBJ_Array);
Map_OBJ_Range=[10 35];
Toral_Intensity_Array=zeros((Map_OBJ_Range(2)-Map_OBJ_Range(1)+1),1);
for p=1:(Map_OBJ_Range(2)-Map_OBJ_Range(1)+1)
    for q=1:length(ROI_LastOBJ_Array)
        if ROI_LastOBJ_Array(q)==(Map_OBJ_Range(1)+p-1)
            Toral_Intensity_Array(p)=Toral_Intensity_Array(p)+ROI_Intensity_at_DET_Array(q);
        end
    end
end
plot(Map_OBJ_Range(1):Map_OBJ_Range(2),Toral_Intensity_Array);
%%
plot(NumberofSegmentforThisRay)