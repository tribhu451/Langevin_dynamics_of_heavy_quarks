import os
import shutil
import subprocess
import re

def generate_ic_input_file(paramsize, paramname, paramval, path) :
  mfile = open("input_parameters", "w")
  mfile.write(
  "dt 0.02  \n"
  "\n\n\n\n"
  "#================================\n"
  )
  mfile.write( "bulk_file %s\n" %(path+"/evolution_xyeta.dat") )
  mfile.write( "ncoll_prof_file %s\n" %(path+"/mc_glauber_boost_invariant_event_averaged_profile_for_rapidity_extension.dat") )
  mfile.write( "glauber_input_file %s\n" %(path+"/init_input") )
  mfile.write( "pythia_file %s\n" %("pythia8309/examples/init_ccbar_info.dat") )
  mfile.close()

  mfile = open("input_parameters", "a")
  for ii in range (0,paramsize) :
    mfile.write( "%s %s\n" %(paramname[ii], paramval[ii]) )
  mfile.close()

  mfile = open("input_parameters", "a")
  mfile.write("\n\nend 1\n\n")
  mfile.close()






#====================================================================================
#====================================================================================
#=================                      MAIN                 ========================
#====================================================================================
#====================================================================================


subprocess.run("ls")
#subprocess.call("rm -rf RUN*", shell = True)
#os.system("rm -rf RUN*")







#============================
#          SET-0            #
#============================

FILENAME = "RUN0"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.420300, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)


#============================
#          SET-1            #
#============================

FILENAME = "RUN1"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.420301, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)


#============================
#          SET-2            #
#============================

FILENAME = "RUN2"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.420302, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)


#============================
#          SET-3            #
#============================

FILENAME = "RUN3"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.420303, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)



#============================
#          SET-4            #
#============================

FILENAME = "RUN4"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.420304, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)



#============================
#          SET-5            #
#============================

FILENAME = "RUN5"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_5060/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.450600, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)


#============================
#          SET-6            #
#============================

FILENAME = "RUN6"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_5060/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.450601, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)


#============================
#          SET-7            #
#============================

FILENAME = "RUN7"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_5060/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.450602, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)



#============================
#          SET-8            #
#============================

FILENAME = "RUN8"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_5060/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.450603, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)



#============================
#          SET-9            #
#============================

FILENAME = "RUN9"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.4_cent_5060/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.450604, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)



#============================
#          SET-10            #
#============================

FILENAME = "RUN10"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.2_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.220300, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)


#============================
#          SET-11            #
#============================

FILENAME = "RUN11"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.2_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.220301, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)



#============================
#          SET-12            #
#============================

FILENAME = "RUN12"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.2_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.220302, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)



#============================
#          SET-13            #
#============================

FILENAME = "RUN13"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.2_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.220303, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)


#============================
#          SET-14            #
#============================

FILENAME = "RUN14"
os.mkdir(FILENAME)
print(FILENAME + " created ... ")
PATH = "/media/tribhuban/Drive3/TRIBHUBAN/HeavyFlavor/DDbar_correlation/LANGEVIN_DYNAMICS_OF_HEAVY_FLAVOR/input_bulk_pbpb_5020gev_tau0_0.2_cent_2030/"
InitparN = [ "Do_evolve_primordial_ccbar_pair", "Do_set_etas_assuming_boost_invariance","delta_etas_gap", "nquarks", "g0", "p0", "t0", "dt" ]
InitparV = [  -1, 1, 0.220304, 20000,  2.0,  0.0, 0.0, 0.02  ]

generate_ic_input_file(len(InitparV), InitparN, InitparV , PATH)
shutil.move("input_parameters", FILENAME)
Aname = ''
for ii in range(0,len(InitparN)) :
    Aname += str( str(InitparN[ii])+ "_" + str(InitparV[ii])+"_" )
ult_name= re.search('lower_(.*)kappa_coefficient', PATH)
xfile = open("file_renaming_and_event_count_info.txt", "w")
supposed_output_filename = Aname
print(supposed_output_filename)
xfile.write(supposed_output_filename + "\n")
xfile.close()
shutil.move("file_renaming_and_event_count_info.txt", FILENAME)











































































































































































