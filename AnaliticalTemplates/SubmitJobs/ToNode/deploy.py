import os
import json
import sys
import shutil
import subprocess

################# Loading JSON dictionary ##################
############################################################

json_file_path = "pars.json"
lPars = list()
currentPath = os.getcwd()

with open(json_file_path, "r") as f_in:
    myPars = json.load(f_in)

#############################################################

######### Creating dir

tmpDir = myPars['tmpJobPath'] + myPars['jobID']
logDir = tmpDir + "/logs"
if not os.path.exists(tmpDir):
    print("\n\n--> Creating tmp Job diretcory...")
    os.mkdir(tmpDir)
    os.mkdir(logDir)
    logDir = logDir + "/"
else:
    print("\n\n--> Directory exists ! Please, check your files and input parameters")
    sys.exit()

######### Building paths

seeds_path = tmpDir + myPars['seeds_file_name']
scaled_path = tmpDir + myPars['scaled_out_path']
templates_allsky = tmpDir + myPars['tmp_templates_allsky']
templates_dampe = tmpDir + myPars['tmp_templates_dampe']

scale_dir = currentPath + myPars['tmp_scale_dir']
template_dir = currentPath + myPars['tmp_template_dir']
seeds_dir = currentPath + myPars['tmp_seeds_dir']
analysis_dir = currentPath + myPars['tmp_analysis_dir']

scale_exe = currentPath + myPars['tmp_scale_exe']
seeds_exe = currentPath + myPars['tmp_seeds_exe']
templates_exe = currentPath + myPars['tmp_templates_exe']

######### Compiling dependencies

print("--> Compiling dependencies...\n\n")

os.chdir(analysis_dir)
bashCommand = "make rebuild"
os.system(bashCommand)

os.chdir(seeds_dir)
bashCommand = "make"
os.system(bashCommand)
bashCommand = "make distclean"
os.system(bashCommand)

os.chdir(scale_dir)
bashCommand = "make"
os.system(bashCommand)
bashCommand = "make distclean"
os.system(bashCommand)

os.chdir(template_dir)
bashCommand = "make rebuild"
os.system(bashCommand)

######### Executing dependencies

print("\n\n--> Executing dependencies...\n\n")

os.chdir(seeds_dir)
temp_command = "{exe:s} {exe_path:s} {ntry:f} {btry:f}"
bashCommand = temp_command.format(exe=seeds_exe,exe_path=seeds_path,ntry=myPars['n_try'],btry=myPars['batch_try'])
os.system(bashCommand)

os.chdir(scale_dir)
temp_command = "{exe:s} {full_reference:s} {scaled_out:s} {Xbinning:f} {Ybinning:f} {LS_events:f} {multi:f}"
bashCommand = temp_command.format(exe=scale_exe,full_reference=myPars['dampe_iso_map'],scaled_out=scaled_path,Xbinning=myPars['binX_scaling'],Ybinning=myPars['binY_scaling'],LS_events=myPars['n_LS'],multi=1)
os.system(bashCommand)

os.chdir(template_dir)
temp_command = "{exe:s} {scaled_maps_path:s} {t_allsky:s} {t_dampe:s} {ldir:s}"
bashCommand = temp_command.format(exe=templates_exe,scaled_maps_path=scaled_path,t_allsky=templates_allsky,t_dampe=templates_dampe,ldir=logDir)
os.system(bashCommand)

######### Launching jobs

print("\n\n--> Launching Jobs...\n\n")

os.chdir(currentPath + "/SHs/")

os.environ["TRY_IDX"] = str(myPars['n_try'])
os.environ["N_TRY"] = str(myPars['n_try'])
os.environ["BATCH_TRY"] = str(myPars['batch_try'])

subprocess.call(['./launcher.sh'])
