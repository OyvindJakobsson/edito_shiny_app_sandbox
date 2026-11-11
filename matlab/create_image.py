import os
import socket
import subprocess

if "dedun" in socket.gethostname():
    # source_files_dirpath = "/home/a23236/algestatus_script_dir/data_root/data_regelmessig/GSM/8-day"
    # images_dirpath = "/home/a23236/algestatus_script_dir/data_root/figs/8-day"
    matlab_root = '/data/osea/scratch/oyvjak/EDITO/matlab'
    # figurer_regelmessig_dirpath = "/home/a23236/algestatus_script_dir/data_root/figurer_regelmessig/"
    matlab_version = "matlab-R2024a"

matlab_script = f"{matlab_version} -nodisplay -nodesktop -nosplash -sd {matlab_root} -batch \"" \
                f"TOPAZ_NRT_CTD_Comp;exit\""

print(matlab_script)
subprocess.run(matlab_script, shell=True)