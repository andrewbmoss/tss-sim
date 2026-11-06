import subprocess,os
with open(r"C:\Users\Andrew\AppData\Local\Temp\imc_probe_a.log","w") as fh:
    r=subprocess.run([r"C:\Users\Andrew\radioconda\envs\tss\python.exe",r"C:\Users\Andrew\AppData\Local\Temp\tss_work\scripts\isis_model_check.py"]+['analyse', 'C:\\Users\\Andrew\\AppData\\Local\\Temp\\imc_probe'],cwd=r"C:\Users\Andrew\AppData\Local\Temp\tss_work",stdout=fh,stderr=subprocess.STDOUT)
    fh.write("\nEXIT %d\n"%r.returncode)
open(r"C:\Users\Andrew\AppData\Local\Temp\imc_probe_a.log.done","w").write("1")
