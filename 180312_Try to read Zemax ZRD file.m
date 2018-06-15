clear all
%%
File_Path='C:\TuanShu\180312_MLOPTIC\180312_SL Test_All Lens 99.ZRD';
% 1st int32: version
% 2nd int32: number of segment (max)
% 3rd in32: number of segment of ray 1
Known_Number_of_Ray=100000;
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
        %Status         fseek(fin,, 'c');
%         unsigned int status;
        RayData{SegNumber,RayNumber}.status=fread(fin,1,'uint32');

%         int level;
        RayData{SegNumber,RayNumber}.level=fread(fin,1,'int32');

%         int hit_object;
        RayData{SegNumber,RayNumber}.hit_object=fread(fin,1,'int32');

%         int hit_face;
        RayData{SegNumber,RayNumber}.hit_face=fread(fin,1,'int32');

%         int unused;
        RayData{SegNumber,RayNumber}.unused=fread(fin,1,'int32');

%         int in_object;
        RayData{SegNumber,RayNumber}.in_object=fread(fin,1,'int32');

%         int parent;
        RayData{SegNumber,RayNumber}.parent=fread(fin,1,'int32');

%         int storage;
        RayData{SegNumber,RayNumber}.storage=fread(fin,1,'int32');

%         int xybin, lmbin;
        RayData{SegNumber,RayNumber}.xybin=fread(fin,1,'int32');
        RayData{SegNumber,RayNumber}.lmbin=fread(fin,1,'int32');

%         double index, starting_phase;
        RayData{SegNumber,RayNumber}.index=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.starting_phase=fread(fin,1,'double');

%         double x, y, z;

        RayData{SegNumber,RayNumber}.x=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.y=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.z=fread(fin,1,'double');

%         double l, m, n;

        RayData{SegNumber,RayNumber}.l=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.m=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.n=fread(fin,1,'double');

%         double nx, ny, nz;

        RayData{SegNumber,RayNumber}.nx=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.ny=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.nz=fread(fin,1,'double');

%         double path_to, intensity;
        RayData{SegNumber,RayNumber}.path_to=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.intensity=fread(fin,1,'double');

%         double phase_of, phase_at;

        RayData{SegNumber,RayNumber}.phase_of=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.phase_at=fread(fin,1,'double');

%         double exr, exi, eyr, eyi, ezr, ezi;

        RayData{SegNumber,RayNumber}.exr=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.exi=fread(fin,1,'double');

        RayData{SegNumber,RayNumber}.eyr=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.eyi=fread(fin,1,'double');

        RayData{SegNumber,RayNumber}.ezr=fread(fin,1,'double');
        RayData{SegNumber,RayNumber}.ezi=fread(fin,1,'double');

        %fseek(fin,Size_of_Segment, 'cof');
        %disp(SegNumber);
    end
    RayNumber=RayNumber+1;
end

%%
SegNumber=29;
RayNumber=9;
RayData{SegNumber,RayNumber}.level