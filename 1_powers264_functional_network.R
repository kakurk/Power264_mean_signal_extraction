#Author: Arnab Roy, Ph.D. Binghamton University, NY-13902
#Post-Doctoral Researcher, Penn. State University, State College-16801
#Date: 24 July 2015

#----------------------------------------------------------------------------------------------------------------------------------

#Objective: The objective of this code is to make connectivity network using powers 264 ROI

#*****************************************************************************************************************************


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

corr_matrix <- function(sig_file,op_name,p_threshold) #corr_matrix.begins
{

library("lsr")

 cormat <- correlate(sig_file,corr.method="pearson",p.adjust.method="fdr")
 cormat$correlation[is.na(cormat$correlation)] <- 0 #for NA cases set correlation to 0
 cormat$p.value[is.na(cormat$p.value)] <- 10000 #For NA cases set the p-value to very large number

 #The NA cases will anyways get excluded as p-value will be large, we will only choose for p-value < 0.05
  
 #-----------------------------------------------------------------------------------

 #Save only the edges that qualify the p-threshold
            
 op_table <- vector('list',(nrow(cormat$correlation)*nrow(cormat$correlation)+1))
 op_table <- list(c('vox_col_id_A','vox_col_id_B','weight','p_value')) #this is the file header

 counter <- 1 #initialize the counter to 1

 for(c4 in 1:(nrow(cormat$correlation)-1)) #For.c4.begins
    {
     for(c5 in (c4+1):ncol(cormat$correlation)) #For.c5.begins
        {
         if(cormat$p.value[c4,c5] < p_threshold)
           {
            vox_col_id_a1 <- c4
            vox_col_id_a2 <- c5

            op_table[[counter+1]] <- c(vox_col_id_a1,vox_col_id_a2,cormat$correlation[c4,c5],cormat$p.value[c4,c5]) 

            counter <- counter + 1

            }

        } #For.c5.ends

    } #For.c4.ends



 #--------------------------------------------------------------------------------------

 file.create(op_name)

 for(filecounter in 1:counter)
    {
     write(c(unlist(op_table[[filecounter]])),file=op_name,append=TRUE,ncol=100)
    }

} #corr_matrix.begins

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


F_powers264_functional_network <- function(IP_f_path,IP_s_path,IP_subject_folder,IP_c_file,IP_bm_file,IP_powers_coordinate,IP_output_file_tag,IP_p_val_threshold)#the main funtion begins
{

#Load the fmri library

library("fmri")

#----------------------------------------------------------------------------------------------------------------

#c1 controls the subjects
for(c1 in 1:length(IP_subject_folder)) #for.c1.begins
   {
    #the functional nifti file
    f_nifti_file_name <- paste(IP_f_path,IP_subject_folder[c1],"/",IP_c_file,sep="")
    print(f_nifti_file_name)
    x <- read.NIFTI(f_nifti_file_name, level = 0.75,setmask=FALSE)
    BOLD_data <- extract.data(x, what = "data")

    #-----------------------------------------------------------------------------------------------------------

    #the structural nifti file
    s_nifti_file_name <- paste(IP_s_path,IP_subject_folder[c1],"/",IP_bm_file[c1],sep="")
    print(s_nifti_file_name)
    mask <- read.NIFTI(s_nifti_file_name, level = 0.75,setmask=FALSE)
    mask_data <- extract.data(mask, what = "data")

    #-----------------------------------------------------------------------------------------------------------

    #if mask dimensions do not match the data dimension exit

    if(x$dim[1] != mask$dim[1] && x$dim[2]  != mask$dim[2] && x$dim[3]  != mask$dim[3])
      {system("echo \"mask dimensions not equal to data dimensions\" >> error ");quit();}


    print(c(x$dim[1],mask$dim[1],x$dim[2],mask$dim[2],x$dim[3],mask$dim[3]))
       
    #-----------------------------------------------------------------------------------------------------------

    #first read powers coordinate file
 
    p_cood <- read.table(IP_powers_coordinate,header=TRUE)
  
    #define empty array for storing the ROI signals
    signal_data <- vector('list',264)  #PARA

    #now iterate through all 264 powers coordinates and find out which voxels belong to each coordinate.
    for(c2 in 1:264)#for.c2.begins  #PARA
       {

        #find the voxel closest to this power's coordinate
        vox_x_cood <- round((91 - p_cood$x[c2])/3)
        vox_y_cood <- round((107 + p_cood$y[c2])/3) 
        vox_z_cood <- round((89 + p_cood$z[c2])/3) 

        #print(c(p_cood$x[c2],p_cood$y[c2],p_cood$z[c2],paste(p_cood$system[c2])))
        #print(c(p_cood$x[c2],p_cood$y[c2],p_cood$z[c2],vox_x_cood,vox_y_cood,vox_z_cood))

 
        #now check all 27 voxels around it and see which ones belong to grey matter.
        #and then create an averge signal using these voxels

        voxel_count <- 0 
        signal <- rep(0,145)  #PARA
        for(c3x in -2:2)#for.c3x.begins
           {
            for(c3y in -2:2)#for.c3y.begins
               {
                for(c3z in -2:2)#for.c3z.begins
                   {

                     if(mask_data[vox_x_cood+c3x,vox_y_cood+c3y,vox_z_cood+c3z,1] > 153)
                       {
                        voxel_count <- voxel_count + 1
                        signal <- signal + BOLD_data[vox_x_cood+c3x,vox_y_cood+c3y,vox_z_cood+c3z,1:145]
                       }

                   }#for.c3z.ends
               }#for.c3y.ends
           }#for.c3x.ends

        #print(c(c2,voxel_count))

        
        #store the average signal
        if(voxel_count == 0)
          {signal_data[[c2]] <- rep(0,145)}  #PARA
        else
          {signal_data[[c2]] <- signal/voxel_count}
      
        
       }#for.c2.ends

    # write out the mean signal for all 264 power ROIs
    outputfilename <- paste(IP_f_path, IP_subject_folder[c1], "/", IP_subject_folder, '_extracted_mean_signal.txt',  sep="")
    lapply(signal_data, write, file=outputfilename, append=TRUE, ncolumns=145)

    #-----------------------------------------------------------------------------------------------------------

    #now create the correlation matrix
    #please note that the signal for each ROI is store across columns.
    #that is col-1 is time-1, col-2 is time-2, and so on.
    #for correlation to work properly we need to transpose this matrix

    signal_data_transposed <- t(do.call(rbind,signal_data))

    corr_matrix(signal_data_transposed,paste(IP_output_file_tag,"_",IP_subject_folder[c1],sep=""),IP_p_val_threshold)
    print("Correlation matrix created.")



   }#for.c1.ends


#----------------------------------------------------------------------------------------------------------------

}#the main funtion ends

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------



functional_path <- '/storage/home/epg5130/scratch/metaconnect/TBI/rest_func/'

structural_path <- '/storage/home/epg5130/scratch/metaconnect/TBI/T1/'

subject_folder <- c('L_britne_hmc_1004', 'L_britne_hmc_1008', 'L_britne_hmc_1009', 'L_britne_hmc_1011', 'L_britne_hmc_1012', 'L_britne_up_1013', 'L_britne_up_1014', 'L_britne_up_1015', 'L_britne_up_1017', 'Long_hmc_1021_1', 'Long_hmc_1022_1', 'Long_hmc_1024_1', 'Long_hmc_1026_1', 'Long_hmc_1027_1', 'Long_hmc_1029_1', 'Long_hmc_1030_1', 'Long_hmc_1031_1', 'Long_hmc_1032_1', 'Long_hmc_1033_1', 'Long_hmc_1036_1', 'Long_hmc_1038_1', 'Long_hmc_1039_1')


conn_output_file<- 'conn_project1_new/results/preprocessing/niftiDATA_Subject001_Condition000.nii'  #PARA

binary_mask_file <- c('resample_wc1cr_co20110630_103138T1MPRAGEIsos013a1001.nii', 'resample_wc1cr_co20110804_161209T1MPRAGEIsos002a1001.nii', 'resample_wc1cr_co20110825_162248T1MPRAGEIsos002a1001.nii', 'resample_wc1cr_co20111117_162457T1MPRAGEIsos002a1001.nii', 'resample_wc1cr_co20111220_164919T1MPRAGEIsos002a1001.nii', 'resample_wc1cr_coT1MPRageSagittal.nii', 'resample_wc1cr_coT1MPRageSagittal.nii', 'resample_wc1cr_coT1MPRageSagittal.nii', 'resample_wc1cr_coT1MPRageSagittals002a1001.nii', 'resample_wc1cr_co20100603_175443T1MPRAGEIsos007a1001.nii', 'resample_wc1cr_co20100610_170305T1MPRAGEIsos007a1001.nii', 'resample_wc1cr_cos012a1001.nii', 'resample_wc1cr_coT1MPRAGEIsos002a1001.nii', 'resample_wc1cr_co20120222_180702T1MPRAGEIsos003a1001A.nii', 'resample_wc1cr_cos002a1001.nii', 'resample_wc1cr_co20120620_171657T1MPRAGEIsos002a1001A.nii', 'resample_wc1cr_co20120718_174845T1MPRAGEIsos002a1001A.nii', 'resample_wc1cr_co20121016_161545T1MPRAGEIsos002a1001.nii', 'resample_wc1cr_co20121016_161545T1MPRAGEIsos002a1001.nii', 'resample_wc1cr_coT1MPRAGEIsos002a1001.nii', 'resample_wc1cr_coT1MPRageSagittal.nii', 'resample_wc1cr_coT1MPRageSagittal.nii')



powers_coordinate <- './powers_MNI_coordinate.dat'  #PARA

output_file_tag <- 'powers_264_conmat_TBI1'  #PARA

p_val_threshold <- 0.05 #PARA


F_powers264_functional_network(functional_path,structural_path,subject_folder,conn_output_file,binary_mask_file,powers_coordinate,output_file_tag,p_val_threshold)






